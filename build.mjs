// @ts-check
import * as esbuild from 'esbuild';

esbuild.buildSync({
  entryPoints: ['src/extension.ts'],
  bundle: true,
  platform: 'node',
  outfile: 'dist/extension.js',
  external: ['vscode', 'node'],
  format: 'cjs',
  minify: true
});
