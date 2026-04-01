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
    return {
        uri: { toString: () => uri, fsPath: '/test.smk' },
        getText: () => text,
        version,
        languageId: 'snakemake',
    };
}

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
// Version cache: update() call count
// ---------------------------------------------------------------------------

test('update() called once for repeated getSymbols() with the same document version', () => {
    const manager = new TrackedManager();
    const doc = makeDoc(SAMPLE, 1);

    manager.getSymbols(doc as any);
    manager.getSymbols(doc as any);
    manager.getSymbols(doc as any);

    assert.equal(manager.updateCount, 1);
});

test('update() called again when document version changes', () => {
    const manager = new TrackedManager();

    manager.getSymbols(makeDoc(SAMPLE, 1) as any);
    manager.getSymbols(makeDoc(SAMPLE, 1) as any);
    manager.getSymbols(makeDoc(SAMPLE, 2) as any);

    assert.equal(manager.updateCount, 2);
});

test('update() called once even across getSymbols(), getLinks(), and checkBlockAt()', () => {
    const manager = new TrackedManager();
    const doc = makeDoc(SAMPLE, 1);

    manager.getSymbols(doc as any);
    manager.getLinks(doc as any);
    // position that doesn't hit any shell block
    manager.checkBlockAt(doc as any, { line: 0, character: 0 } as any);

    assert.equal(manager.updateCount, 1);
});

test('independent documents have separate caches', () => {
    const manager = new TrackedManager();
    const docA = makeDoc(SAMPLE, 1, 'file:///a.smk');
    const docB = makeDoc(SAMPLE, 1, 'file:///b.smk');

    manager.getSymbols(docA as any);
    manager.getSymbols(docA as any); // cache hit for A
    manager.getSymbols(docB as any); // cache miss for B
    manager.getSymbols(docB as any); // cache hit for B

    assert.equal(manager.updateCount, 2); // one per document
});

// ---------------------------------------------------------------------------
// Version cache: correctness of returned data
// ---------------------------------------------------------------------------

test('returns the same symbols array reference on cache hit', () => {
    const manager = new SnakemakeManager();
    const doc = makeDoc(SAMPLE, 1);

    const first = manager.getSymbols(doc as any);
    const second = manager.getSymbols(doc as any);

    assert.strictEqual(first, second);
});

test('correctly extracts rule names from snakemake document', () => {
    const manager = new SnakemakeManager();
    const symbols = manager.getSymbols(makeDoc(SAMPLE, 1) as any);

    assert.equal(symbols.length, 2);
    assert.equal(symbols[0].name, 'foo');
    assert.equal(symbols[1].name, 'bar');
});

test('symbols update when document version changes', () => {
    const manager = new SnakemakeManager();

    const v1 = manager.getSymbols(makeDoc(SAMPLE, 1) as any);
    assert.equal(v1.length, 2);

    const extendedSample = SAMPLE + '\nrule baz:\n    input: "c.txt"\n';
    const v2 = manager.getSymbols(makeDoc(extendedSample, 2) as any);
    assert.equal(v2.length, 3);
    assert.equal(v2[2].name, 'baz');
});

// ---------------------------------------------------------------------------
// Debounce: skips redundant update when version already parsed
// ---------------------------------------------------------------------------

test('debounced update skips if getSymbols() already parsed that version', async () => {
    const manager = new TrackedManager();
    const doc = makeDoc(SAMPLE, 1);

    // Simulate provider call that parses the document synchronously.
    manager.getSymbols(doc as any);
    assert.equal(manager.updateCount, 1);

    // Trigger the debounce (as VS Code would on document open/change).
    manager.onDocumentChange(doc as any);

    // Wait for the 200 ms debounce to fire.
    await new Promise(resolve => setTimeout(resolve, 250));

    // The debounce callback should have seen versionMap already up-to-date and skipped.
    assert.equal(manager.updateCount, 1);
});

test('debounced update fires when version has not yet been parsed', async () => {
    const manager = new TrackedManager();
    const doc1 = makeDoc(SAMPLE, 1);
    const doc2 = makeDoc(SAMPLE + '\nrule baz:\n    input: "c.txt"\n', 2, 'file:///test.smk');

    // Parse v1.
    manager.getSymbols(doc1 as any);
    assert.equal(manager.updateCount, 1);

    // Document changes to v2 – trigger the debounce but do NOT call getSymbols.
    manager.onDocumentChange(doc2 as any);

    await new Promise(resolve => setTimeout(resolve, 250));

    // Debounce must have fired because v2 had not yet been parsed.
    assert.equal(manager.updateCount, 2);
});

test('rapid onDocumentChange calls result in only one debounced update', async () => {
    const manager = new TrackedManager();

    // Simulate initial parse so we have a version baseline.
    manager.getSymbols(makeDoc(SAMPLE, 1) as any);
    assert.equal(manager.updateCount, 1);

    // Fire onDocumentChange many times in quick succession (v2 doc).
    const doc2 = makeDoc(SAMPLE, 2);
    for (let i = 0; i < 10; i++) {
        manager.onDocumentChange(doc2 as any);
    }

    await new Promise(resolve => setTimeout(resolve, 250));

    // Only one debounced update should have fired.
    assert.equal(manager.updateCount, 2);
});
