import * as vscode from 'vscode';
import * as path from 'path';
import * as fs from 'fs';


export function activate(context: vscode.ExtensionContext) {
    const snakemakeManager = new SnakemakeManager();
    const selector = { scheme: 'file', language: 'snakemake' }

    context.subscriptions.push(
        vscode.workspace.onDidChangeTextDocument(snakemakeManager.onDocumentChange),
        vscode.languages.registerHoverProvider(selector, {
            provideHover(document, position, token) {
                const block = snakemakeManager.checkBlockAt(document, position);
                if (!block) return;

                const md = new vscode.MarkdownString();
                md.appendCodeblock(trimEmptyLines(block.content), 'shellscript');

                return new vscode.Hover(md, block.range);
            }
        }),
        vscode.languages.registerDocumentSymbolProvider(selector, {
            provideDocumentSymbols(document: vscode.TextDocument) {
                return snakemakeManager.getSymbols(document);
            }
        }),
        vscode.languages.registerDocumentLinkProvider(selector, {
            provideDocumentLinks(document: vscode.TextDocument) {
                return snakemakeManager.getLinks(document);
            }
        })
    );
}


interface ShellBlock {
    range: vscode.Range;
    content: string;
}

const snakemakeRegex = {
    rule: /^\s*(rule|checkpoint)(?:\s+([a-zA-Z_][a-zA-Z0-9_]*))?\s*:/,
    module: /^\s*(module|subworkflow)(?:\s+([a-zA-Z_][a-zA-Z0-9_]*))?\s*:/,
    shell: /^\s+shell\s*:/,
    path: {
        global: /^\s*(conda|include|configfile)\s*:/,
        rule: /^\s+(conda|script)\s*:/,
        module: /^\s+(snakefile)\s*:/,
    }
}

class SnakemakeManager {
    private blockMap = new Map<string, ShellBlock[]>();
    private symbolMap = new Map<string, vscode.DocumentSymbol[]>();
    private linkMap = new Map<string, vscode.DocumentLink[]>();
    private updateTimeout = new Map<string, NodeJS.Timeout>();

    onDocumentChange(event: vscode.TextDocumentChangeEvent) {
        if (event.document.languageId.toLowerCase() !== 'snakemake') return;

        if (this.updateTimeout.has(event.document.uri.toString())) {
            clearTimeout(this.updateTimeout.get(event.document.uri.toString())!);
        }
        this.updateTimeout.set(event.document.uri.toString(), setTimeout(() => {
            this.update(event.document);
            this.updateTimeout.delete(event.document.uri.toString());
        }, 200));
    }

    update(document: vscode.TextDocument) {
        const blocks: ShellBlock[] = [];
        const symbols: vscode.DocumentSymbol[] = [];
        const links: vscode.DocumentLink[] = [];
        this.blockMap.set(document.uri.toString(), blocks);
        this.symbolMap.set(document.uri.toString(), symbols);
        this.linkMap.set(document.uri.toString(), links);

        const text = document.getText();
        const lines = text.split('\n');

        let state = {
            shellBlock: { delimiter: "", content: [''] },
            blockStart: new vscode.Position(0, 0),
            current: { indent: 0, blockType: 'global' },
        };
        state.shellBlock.content = [];

        let regexMatch: RegExpExecArray | null;

        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];
            if (state.shellBlock.delimiter) {
                this.pushAnyBlocks(blocks, state, i, line, 0);
                continue;
            }

            if (!line.trim() || line.trim().startsWith('#')) continue;

            const indentMatch = line.match(/^(\s*)[\S]/);
            const indent = indentMatch ? indentMatch[1].length : 0;
            if (indent <= state.current.indent) {
                state.current = { indent: indent, blockType: 'global' };
            }
            if (regexMatch = /(.+):path$/.exec(state.current.blockType)) {
                const link = this.checkLinks(line, i, document.uri.fsPath)
                if (link) {
                    links.push(link)
                }
                // the attempt failed, never try anymore
                state.current.blockType = regexMatch[1]
                continue
            }
            if (state.current.blockType === 'global') {
                if (regexMatch = snakemakeRegex.rule.exec(line)) {
                    const ruleType = regexMatch[1];
                    const ruleName = regexMatch[2] || 'unnamed';
                    const range = new vscode.Range(i, 0, i, line.length);
                    const symbol = new vscode.DocumentSymbol(ruleName, ruleType, vscode.SymbolKind.Function, range, range);
                    symbols.push(symbol);
                    state.current.blockType = 'rule'
                    continue;
                } else if (regexMatch = snakemakeRegex.module.exec(line)) {
                    const ruleType = regexMatch[1];
                    const ruleName = regexMatch[2] || 'unnamed';
                    const range = new vscode.Range(i, 0, i, line.length);
                    const symbol = new vscode.DocumentSymbol(ruleName, ruleType, vscode.SymbolKind.Function, range, range);
                    symbols.push(symbol);
                    state.current.blockType = 'module'
                    continue;
                } else if (
                    snakemakeRegex.path.global.exec(line)
                ) {
                    state.current.blockType = 'global:path'
                    // } else if (
                    //     snakemakeRegex.shell.test(line) ||
                    //     snakemakeRegex.path.rule.test(line) ||
                    //     snakemakeRegex.path.module.test(line)
                    // ) {
                    //     error
                }
            } else if (state.current.blockType === 'module') {
                if (snakemakeRegex.path.module.exec(line)) {
                    state.current.blockType = 'module:path'
                }
            } else if (state.current.blockType === 'rule') {
                if (snakemakeRegex.shell.test(line)) {
                    state.current.blockType = 'shell';
                } else if (snakemakeRegex.path.rule.exec(line)) {
                    state.current.blockType = 'rule:path'
                }
            }
            if (regexMatch = /(.+):path$/.exec(state.current.blockType)) {
                const link = this.checkLinks(line, i, document.uri.fsPath)
                if (link) {
                    links.push(link)
                    // check the file path just after the keyword
                    state.current.blockType = regexMatch[1]
                }
            } else if (state.current.blockType === 'shell') {
                const shellMatch = line.match(/("""|''')/);
                if (shellMatch) {
                    state.shellBlock = { delimiter: shellMatch[0], content: [] };
                    state.blockStart = new vscode.Position(i, line.indexOf(state.shellBlock.delimiter) + state.shellBlock.delimiter.length);
                }
                this.pushAnyBlocks(blocks, state, i, line, state.blockStart.character);
            }
        }
    }

    private pushAnyBlocks(
        blocks: ShellBlock[],
        state = {
            shellBlock: { delimiter: "", content: [''] },
            blockStart: new vscode.Position(0, 0)
        },
        i = 0, line = '', shift = 0
    ) {
        while (state.shellBlock.delimiter) {
            const endDelimiterIndex = line.indexOf(state.shellBlock.delimiter, shift);

            if (endDelimiterIndex === -1) {
                state.shellBlock.content.push(line.substring(shift));
                break;
            }
            state.shellBlock.content.push(line.substring(shift, shift + endDelimiterIndex));


            const blockRange = new vscode.Range(state.blockStart, new vscode.Position(i, endDelimiterIndex));
            blocks.push({ range: blockRange, content: state.shellBlock.content.join('\n') });

            const newStart = endDelimiterIndex + state.shellBlock.delimiter.length;
            const newStartMatch = line.substring(newStart).match(/("""|''')/);
            if (newStartMatch) {
                state.shellBlock = { delimiter: newStartMatch[0], content: [] };
                const newDelimiterPos = line.indexOf(newStartMatch[0], newStart);
                state.blockStart = new vscode.Position(i, newDelimiterPos + newStartMatch[0].length);
            } else {
                state.shellBlock = { delimiter: "", content: [] };
            }
        }
    }

    checkLinks(line: string, i: number, fsPath: string) {
        if (/^[^'"]*#/.test(line)) {
            return
        }
        const match = /\s*("([^"]+)")|('([^']+)')/.exec(line);
        if (match) {
            const filePath = match[2] || match[4];
            const startCol = line.indexOf(filePath);
            const endCol = startCol + filePath.length;

            const baseDir = path.dirname(fsPath);
            const fullPath = path.resolve(baseDir, filePath);

            if (fs.existsSync(fullPath)) {
                const linkRange = new vscode.Range(
                    new vscode.Position(i, startCol),
                    new vscode.Position(i, endCol)
                );
                const targetUri = vscode.Uri.file(fullPath);
                return new vscode.DocumentLink(linkRange, targetUri);
            }
        }
    }

    checkDocument(document: vscode.TextDocument) {
        const documentNeedsUpdate = !this.blockMap.get(document.uri.toString());
        if (documentNeedsUpdate) {
            this.update(document);
        }
    }

    checkBlockAt(document: vscode.TextDocument, position: vscode.Position): ShellBlock | undefined {
        this.checkDocument(document);
        const blocks = this.blockMap.get(document.uri.toString()) ?? [];
        return blocks.find(b => b.range.contains(position));
    }

    getSymbols(document: vscode.TextDocument): vscode.DocumentSymbol[] {
        this.checkDocument(document);
        return this.symbolMap.get(document.uri.toString()) || [];
    }

    getLinks(document: vscode.TextDocument): vscode.DocumentLink[] {
        this.checkDocument(document);
        return this.linkMap.get(document.uri.toString()) || [];
    }

}


function trimEmptyLines(block: string): string {
    const lines = block.split('\n');

    let start = 0;
    while (start < lines.length && lines[start].trim() === '') {
        start++;
    }

    let end = lines.length - 1;
    while (end >= 0 && lines[end].trim() === '') {
        end--;
    }

    return lines.slice(start, end + 1).join('\n');
}
