// @ts-check
import * as esbuild from 'esbuild';
import { resolve, dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));

esbuild.buildSync({
    entryPoints: ['src/test/snakemake-manager.test.ts'],
    bundle: true,
    platform: 'node',
    outfile: 'dist/test/snakemake-manager.test.js',
    // Replace the real 'vscode' module with our lightweight mock.
    alias: { vscode: resolve(__dirname, 'src/test/vscode-mock.ts') },
    format: 'cjs',
    // Keep Node built-ins and the test runner external (resolved at runtime).
    external: ['node:test', 'node:assert', 'node:assert/strict', 'fs', 'path'],
});
