import { test } from 'node:test';
import assert from 'node:assert/strict';
import { SnakemakeManager } from '../extension';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

interface MockDocument {
    uri: { toString(): string; fsPath: string };
    getText(): string;
    version: number;
    languageId: string;
}

function makeDoc(text: string, version: number, uri = 'file:///test.smk'): MockDocument {
    const fsPath = uri.replace(/^file:\/\//, '');
    return {
        uri: { toString: () => uri, fsPath },
        getText: () => text,
        version,
        languageId: 'snakemake',
    };
}

const POS = { line: 0, character: 0 };

/** Subclass that counts how many times update() is actually invoked. */
class TrackedManager extends SnakemakeManager {
    updateCount = 0;
    override update(document: Parameters<SnakemakeManager['update']>[0]) {
        this.updateCount++;
        super.update(document);
    }
}

const SAMPLE = `
rule foo:
    input: "a.txt"
    output: "b.txt"
    shell: """
        echo hello
    """

rule bar:
    input: "b.txt"
    shell: """
        echo world
    """
`;

// ---------------------------------------------------------------------------
// Core invariant: getSymbols / getLinks never trigger update()
// ---------------------------------------------------------------------------

test('getSymbols() never triggers update()', () => {
    const manager = new TrackedManager();
    const doc = makeDoc(SAMPLE, 1);

    manager.getSymbols(doc as any);
    manager.getSymbols(doc as any);
    manager.getSymbols(doc as any);

    assert.equal(manager.updateCount, 0);
});

test('getLinks() never triggers update()', () => {
    const manager = new TrackedManager();
    manager.getLinks(makeDoc(SAMPLE, 1) as any);
    manager.getLinks(makeDoc(SAMPLE, 1) as any);
    assert.equal(manager.updateCount, 0);
});

// ---------------------------------------------------------------------------
// checkBlockAt: sync parse only on first access per document
// ---------------------------------------------------------------------------

test('checkBlockAt() triggers update() on first access, not on subsequent calls', () => {
    const manager = new TrackedManager();
    const doc = makeDoc(SAMPLE, 1);

    manager.checkBlockAt(doc as any, POS as any);
    manager.checkBlockAt(doc as any, POS as any);
    manager.checkBlockAt(doc as any, POS as any);

    assert.equal(manager.updateCount, 1);
});

test('checkBlockAt() does not re-parse when version changes (stale data is acceptable)', () => {
    const manager = new TrackedManager();
    // First access parses v1.
    manager.checkBlockAt(makeDoc(SAMPLE, 1) as any, POS as any);
    // Version changed, but no re-parse – debounce handles this.
    manager.checkBlockAt(makeDoc(SAMPLE, 2) as any, POS as any);

    assert.equal(manager.updateCount, 1);
});

test('checkBlockAt() parses each new document URI independently', () => {
    const manager = new TrackedManager();
    manager.checkBlockAt(makeDoc(SAMPLE, 1, 'file:///a.smk') as any, POS as any);
    manager.checkBlockAt(makeDoc(SAMPLE, 1, 'file:///b.smk') as any, POS as any);
    assert.equal(manager.updateCount, 2);
});

// ---------------------------------------------------------------------------
// Correctness: symbols/links populated after explicit update()
// ---------------------------------------------------------------------------

test('correctly extracts rule names after update()', () => {
    const manager = new SnakemakeManager();
    const doc = makeDoc(SAMPLE, 1);
    manager.update(doc as any);

    const symbols = manager.getSymbols(doc as any);
    assert.equal(symbols.length, 2);
    assert.equal(symbols[0].name, 'foo');
    assert.equal(symbols[1].name, 'bar');
});

test('getSymbols() returns the same array reference on repeated calls (no re-allocation)', () => {
    const manager = new SnakemakeManager();
    const doc = makeDoc(SAMPLE, 1);
    manager.update(doc as any);

    const first = manager.getSymbols(doc as any);
    const second = manager.getSymbols(doc as any);
    assert.strictEqual(first, second);
});

test('symbols reflect new content after update() called with a new version', () => {
    const manager = new SnakemakeManager();
    const doc1 = makeDoc(SAMPLE, 1);
    manager.update(doc1 as any);
    assert.equal(manager.getSymbols(doc1 as any).length, 2);

    const extendedSample = SAMPLE + '\nrule baz:\n    input: "c.txt"\n';
    const doc2 = makeDoc(extendedSample, 2);
    manager.update(doc2 as any);
    const symbols = manager.getSymbols(doc2 as any);
    assert.equal(symbols.length, 3);
    assert.equal(symbols[2].name, 'baz');
});

// ---------------------------------------------------------------------------
// Debounce: onDocumentChange drives all updates during normal editing
// ---------------------------------------------------------------------------

test('onDocumentChange() populates symbols after debounce fires', async () => {
    const manager = new TrackedManager();
    const doc = makeDoc(SAMPLE, 1);

    assert.equal(manager.getSymbols(doc as any).length, 0); // empty before any parse

    manager.onDocumentChange(doc as any);
    await new Promise(resolve => setTimeout(resolve, 250));

    assert.equal(manager.updateCount, 1);
    assert.equal(manager.getSymbols(doc as any).length, 2);
});

test('debounced update skips when version was already parsed by checkBlockAt()', async () => {
    const manager = new TrackedManager();
    const doc = makeDoc(SAMPLE, 1);

    // Hover forces a sync parse on first access.
    manager.checkBlockAt(doc as any, POS as any);
    assert.equal(manager.updateCount, 1);

    // VS Code also fires onDocumentChange for the same version.
    manager.onDocumentChange(doc as any);
    await new Promise(resolve => setTimeout(resolve, 250));

    // Debounce should see versionMap already up-to-date and skip.
    assert.equal(manager.updateCount, 1);
});

test('debounced update fires when a new version has not yet been parsed', async () => {
    const manager = new TrackedManager();
    const doc1 = makeDoc(SAMPLE, 1);
    const extendedSample = SAMPLE + '\nrule baz:\n    input: "c.txt"\n';
    const doc2 = makeDoc(extendedSample, 2);

    // Parse v1 via debounce.
    manager.onDocumentChange(doc1 as any);
    await new Promise(resolve => setTimeout(resolve, 250));
    assert.equal(manager.updateCount, 1);
    assert.equal(manager.getSymbols(doc1 as any).length, 2);

    // Document changes to v2 – trigger debounce but no provider call.
    manager.onDocumentChange(doc2 as any);
    await new Promise(resolve => setTimeout(resolve, 250));

    assert.equal(manager.updateCount, 2);
    assert.equal(manager.getSymbols(doc2 as any).length, 3);
});

test('rapid onDocumentChange() calls result in only one debounced update', async () => {
    const manager = new TrackedManager();

    // Establish v1 baseline.
    manager.onDocumentChange(makeDoc(SAMPLE, 1) as any);
    await new Promise(resolve => setTimeout(resolve, 250));
    assert.equal(manager.updateCount, 1);

    // Simulate fast typing: 10 changes in quick succession, all version 2.
    const doc2 = makeDoc(SAMPLE, 2);
    for (let i = 0; i < 10; i++) {
        manager.onDocumentChange(doc2 as any);
    }
    await new Promise(resolve => setTimeout(resolve, 250));

    assert.equal(manager.updateCount, 2); // only one additional parse
});
