import { mkdir, readFile, rm, writeFile } from 'node:fs/promises';

import { build } from 'esbuild';

const entryPoint = 'src/polynomial.ts';
const dependencies = ['complex.js', 'fraction.js', 'quaternion'];
const shared = {
  bundle: true,
  entryPoints: [entryPoint],
  legalComments: 'inline',
  target: 'es2020',
};

await mkdir('dist', { recursive: true });

await Promise.all([
  build({
    ...shared,
    external: dependencies,
    footer: {
      js: [
        'const PolynomialExport = module.exports.default;',
        'Object.defineProperty(PolynomialExport, "__esModule", { value: true });',
        'PolynomialExport.default = PolynomialExport;',
        'PolynomialExport.Polynomial = PolynomialExport;',
        'module.exports = PolynomialExport;',
      ].join('\n'),
    },
    format: 'cjs',
    outfile: 'dist/polynomial.js',
    platform: 'node',
    sourcemap: true,
  }),
  build({
    ...shared,
    external: dependencies,
    format: 'esm',
    outfile: 'dist/polynomial.mjs',
    sourcemap: true,
  }),
  build({
    ...shared,
    footer: {
      js: [
        'var PolynomialExport = PolynomialExports.default;',
        'if (typeof define === "function" && define.amd) define([], function () { return PolynomialExport; });',
        'else (typeof globalThis !== "undefined" ? globalThis : typeof self !== "undefined" ? self : this).Polynomial = PolynomialExport;',
      ].join('\n'),
    },
    format: 'iife',
    globalName: 'PolynomialExports',
    minify: true,
    outfile: 'dist/polynomial.min.js',
  }),
  build({
    ...shared,
    format: 'esm',
    minify: true,
    outfile: 'dist/polynomial.min.mjs',
  }),
]);

const declaration = await readFile('.types/polynomial.d.ts', 'utf8');
const commonJsDeclaration = declaration
  .replaceAll('PolynomialResidueInput<', 'PolynomialResidueInputRoot<')
  .replaceAll('PolynomialQuotientRing<', 'PolynomialQuotientRingRoot<')
  .replaceAll('PolynomialResidue<', 'PolynomialResidueRoot<')
  .replaceAll('BuiltInField', 'BuiltInFieldRoot')
  .replaceAll('PolynomialField<', 'PolynomialFieldRoot<')
  .replaceAll('CoefficientMap<', 'CoefficientMapRoot<')
  .replaceAll('PolynomialInput<', 'PolynomialInputRoot<')
  .replace('export type BuiltInFieldRoot', 'type BuiltInFieldRoot')
  .replace('export interface PolynomialFieldRoot', 'interface PolynomialFieldRoot')
  .replace('export type CoefficientMapRoot', 'type CoefficientMapRoot')
  .replace('export type PolynomialInputRoot', 'type PolynomialInputRoot')
  .replace('export type PolynomialResidueInputRoot', 'type PolynomialResidueInputRoot')
  .replace('export interface PolynomialQuotientRingRoot', 'interface PolynomialQuotientRingRoot')
  .replace('export interface PolynomialResidueRoot', 'interface PolynomialResidueRoot')
  .replace('export declare class Polynomial', 'declare class Polynomial')
  .replace(
    'export default Polynomial;',
    `declare namespace Polynomial {
  type FieldName = BuiltInFieldRoot;
  type Field<C> = PolynomialFieldRoot<C>;
  type Coefficients<C> = CoefficientMapRoot<C>;
  type Input<C> = PolynomialInputRoot<C>;
  type QuotientRing<C> = PolynomialQuotientRingRoot<C>;
  type Residue<C> = PolynomialResidueRoot<C>;
  type ResidueInput<C> = PolynomialResidueInputRoot<C>;
}
export = Polynomial;`,
  );

await Promise.all([
  writeFile('dist/polynomial.d.mts', declaration),
  writeFile('dist/polynomial.d.ts', commonJsDeclaration),
]);
await rm('.types', { recursive: true, force: true });