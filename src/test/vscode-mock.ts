/**
 * Minimal mock of the VS Code API for unit testing SnakemakeManager
 * without a running VS Code instance.
 */

export class Position {
    constructor(public readonly line: number, public readonly character: number) {}
}

export class Range {
    readonly start: Position;
    readonly end: Position;

    constructor(startOrLine: Position | number, endOrChar: Position | number, endLine?: number, endChar?: number) {
        if (startOrLine instanceof Position) {
            this.start = startOrLine;
            this.end = endOrChar as Position;
        } else {
            this.start = new Position(startOrLine as number, endOrChar as number);
            this.end = new Position(endLine!, endChar!);
        }
    }

    contains(position: Position): boolean {
        if (position.line < this.start.line || position.line > this.end.line) return false;
        if (position.line === this.start.line && position.character < this.start.character) return false;
        if (position.line === this.end.line && position.character > this.end.character) return false;
        return true;
    }
}

export enum SymbolKind {
    Function = 11,
}

export class DocumentSymbol {
    constructor(
        public name: string,
        public detail: string,
        public kind: SymbolKind,
        public range: Range,
        public selectionRange: Range,
    ) {}
}

export class DocumentLink {
    constructor(public range: Range, public target?: { toString(): string }) {}
}

export const Uri = {
    file(p: string) {
        return { toString: () => `file://${p}`, fsPath: p };
    }
};

export class MarkdownString {
    appendCodeblock(_code: string, _language?: string) {}
}

export class Hover {
    constructor(public contents: MarkdownString, public range?: Range) {}
}
