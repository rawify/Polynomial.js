'use strict';

const Polynomial = require('polynomial');

const bigintField = {
  add: (left, right) => BigInt(left) + BigInt(right),
  sub: (left, right) => BigInt(left) - BigInt(right),
  mul: (left, right) => BigInt(left) * BigInt(right),
  div: (left, right) => BigInt(left) / BigInt(right),
  parse: value => BigInt(value),
  empty: value => value === undefined || BigInt(value) === 0n,
  pow: (base, exponent) => BigInt(base) ** BigInt(exponent),
  equals: (left, right) => BigInt(left) === BigInt(right),
};

const polynomial = new Polynomial('3891238912389182391919282829x^3+4', bigintField);
console.log(polynomial.add('23819218x^2-381').toString());
