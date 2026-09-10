'use strict';

const assert = require('node:assert/strict');
const { readFile } = require('node:fs/promises');
const { describe, it } = require('node:test');
const vm = require('node:vm');

describe('package distributions', function () {
  it('exports a constructable CommonJS class', function () {
    const Polynomial = require('../dist/polynomial.js');

    const polynomial = new Polynomial('x^2-1');
    assert.equal(polynomial.toString(), 'x^2-1');
    assert.equal(Polynomial.default, Polynomial);
    assert.equal(Polynomial.Polynomial, Polynomial);
  });

  it('exports the same class as ESM default and named exports', async function () {
    const module = await import('../dist/polynomial.mjs');

    assert.equal(module.default, module.Polynomial);
    assert.equal(new module.Polynomial('x+1').toString(), 'x+1');
  });

  it('provides a standalone browser global', async function () {
    const source = await readFile('dist/polynomial.min.js', 'utf8');
    const context = vm.createContext({});

    vm.runInContext(source, context);

    assert.equal(typeof context.Polynomial, 'function');
    context.Polynomial.setField('Q');
    assert.equal(new context.Polynomial('1/3x').add('1/3x').toString(), '0.(6)x');
    assert.equal(new context.Polynomial('x^2+1').compose('x+1').toString(), 'x^2+2x+2');
    assert.equal(new context.Polynomial('x^2-3x').dispersion(), 3);
    assert.equal(new context.Polynomial('9x^2-30x+25').sqrt().toString(), '3x-5');
    context.Polynomial.setField('Z2');
    const quotientRing = context.Polynomial.quotient('x^3+x+1');
    assert.equal(quotientRing.create('x').pow(7).toString(), '1');
    assert.equal(quotientRing.isRoot('x^3+x+1', quotientRing.create('x')), true);
    assert.equal(quotientRing.create('x').minimalPolynomial().toString(), 'x^3+x+1');
    assert.equal(quotientRing.primitiveElements().length, 6);
    assert.equal(
      quotientRing.primitivePolynomials().map(String).join(','),
      'x^3+x+1,x^3+x^2+1',
    );
    assert.equal(quotientRing.create('x').multiplicativeSubgroup().length, 7);
    assert.equal(
      new context.Polynomial('x^4+x^2+1').factor().map(String).join(','),
      'x^2+x+1,x^2+x+1',
    );
    assert.equal(
      new context.Polynomial('x^2+x+1').cyclicCodeBasis(3).map(String).join(','),
      'x^2+x+1',
    );
    assert.equal(
      context.Polynomial.cyclicCodeGenerators(7, 4).map(String).join(','),
      'x^3+x^2+1,x^3+x+1',
    );
    assert.equal(context.Polynomial.enumerate(2).length, 8);
  });

  it('provides a standalone AMD module', async function () {
    const source = await readFile('dist/polynomial.min.js', 'utf8');
    let exported;
    const define = (dependencies, factory) => {
      assert.equal(dependencies.length, 0);
      exported = factory();
    };
    define.amd = {};

    vm.runInNewContext(source, { define });

    assert.equal(typeof exported, 'function');
    assert.equal(new exported('x^2').degree(), 2);
  });

  it('provides a standalone browser ESM build', async function () {
    const module = await import('../dist/polynomial.min.mjs');

    module.default.setField('C');
    assert.equal(new module.default('x-i').mul('x+i').toString(), 'x^2+1');
  });
});