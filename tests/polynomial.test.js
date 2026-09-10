'use strict';

const assert = require('node:assert/strict');
const { afterEach, describe, it } = require('node:test');

const Poly = require('../dist/polynomial.js');

function modulo(value, modulus) {
  return ((Number(value) % modulus) + modulus) % modulus;
}

function createModuloField(modulus) {
  return {
    add: (left, right) => modulo(left + Number(right), modulus),
    sub: (left, right) => modulo(left - Number(right), modulus),
    mul: (left, right) => modulo(left * Number(right), modulus),
    div: (left, right) => {
      for (let candidate = 1; candidate < modulus; candidate += 1) {
        if (modulo(Number(right) * candidate, modulus) === 1) {
          return modulo(left * candidate, modulus);
        }
      }
      throw new RangeError('Coefficient is not invertible');
    },
    parse: value => modulo(value, modulus),
    empty: value => value === undefined || value === 0,
    pow: (base, exponent) => modulo(base ** exponent, modulus),
    equals: (left, right) => left === modulo(right, modulus),
  };
}

afterEach(function () {
  Poly.trace = null;
  Poly.setField('R');
});

var tests = [
  {f: "C", o: "add", a: "-1x^3-1x^2-1x-1", b: 0, r: "-x^3-x^2-x-1"},
  {f: "R", o: "add", a: "-1x^3-1x^2-1x-1", b: 0, r: "-x^3-x^2-x-1"},
  {f: "Q", o: "add", a: "-1x^3-1x^2-1x-1", b: 0, r: "-x^3-x^2-x-1"},
  {f: "C", o: "add", a: "+1x^3+1x^2+1x+1", b: 0, r: "x^3+x^2+x+1"},
  {f: "R", o: "add", a: "+1x^3+1x^2+1x+1", b: 0, r: "x^3+x^2+x+1"},
  {f: "Q", o: "add", a: "+1x^3+1x^2+1x+1", b: 0, r: "x^3+x^2+x+1"},
  {f: 'R', o: 'mul', a: undefined, b: null, r: "0"},
  {f: 'R', o: 'mul', a: [-3, 5, -2], b: [-3, 5, 7], r: "-14x^4+25x^3+10x^2-30x+9"},
  {f: 'Z7', o: 'mul', a: [-2, 5], b: [-1, 3, 4], r: "6x^3+3x+2"},
  {f: 'Z7', o: 'div', a: [3, 1, 1, 5, 3, 1], b: [0, 1, 4, 2], r: "4x^2+4x+3"},
  {f: 'Z7', o: 'mod', a: [3, 1, 1, 5, 3, 1], b: [0, 1, 4, 2], r: "6x^2+5x+3"},
  {f: 'Q', o: 'div', a: [3, 1, 1, 5, 3, 1], b: [0, 1, 4, 2], r: "0.5x^2+0.5x+1.25"},
  {f: 'Q', o: 'mod', a: [0, 1, 4, 2], b: [3, 5, 6], r: "-1.9(4)x-1.1(6)"},
  {f: 'R', o: 'add', a: "22-14x^4+25x^3+10x^2-30x+9", b: "-2+1000x", r: "-14x^4+25x^3+10x^2+970x+29"},
  {f: 'Q', o: 'sub', a: "-1.9(4)x-1.1(6)", b: -3, r: "-1.9(4)x+1.8(3)"},
  {f: 'Q', o: 'addmul', a: 2, b: 5, c: 6, r: "32"},
  {f: 'Q', o: 'add', a: "1x", b: 2, r: "x+2"},
  {f: 'Q', o: 'derive', a: "5+3x^3+6x^5", b: 0, r: "6x^5+3x^3+5"},
  {f: 'Q', o: 'derive', a: "5+3x^3+6x^5", b: 1, r: "30x^4+9x^2"},
  {f: 'Q', o: 'derive', a: "5+3x^3+6x^5", b: 2, r: "120x^3+18x"},
  {f: 'Q', o: 'derive', a: "92x+3+5x^2+43x^4", b: 2, r: "516x^2+10"},
  {f: 'Q', o: 'derive', a: "1+2x+3x^2+4x^3+5x^4+6x^5+7x^6", b: 3, r: "840x^3+360x^2+120x+24"},
  {f: 'Q', o: 'derive', a: "1+2x+3x^2+4x^3+5x^4+6x^5+7x^6", b: 0, r: "7x^6+6x^5+5x^4+4x^3+3x^2+2x+1"},
  {f: 'Z2', o: 'sub', a: 2, b: 1, r: '1'},
  {f: "C", o: "mul", a: "x-i", b: "x+i", r: "x^2+1"},
  {f: 'R', o: 'add', a: 0, b: 0, r: '0'},
  {f: 'R', o: 'add', a: 2, b: 3, r: '5'},
  {f: 'Q', o: 'add', a: '0.(3)x', b: '0.(6)x', r: 'x'},
  {f: "Q", o: "div", a: "x^5+3x^4+5x^3+x^2+x+3", b: "2x^3+4x^2+x", r: "0.5x^2+0.5x+1.25"},
  {f: "Q", o: "mod", a: "x^5+3x^4+5x^3+x^2+x+3", b: "2x^3+4x^2+x", r: "-4.5x^2-0.25x+3"},
  {f: 'R', o: 'sub', a: 2, b: 3, r: '-1'},
  {f: 'R', o: 'div', a: "x^4-81", b: "3x^2+27", r: '0.3333333333333333x^2-3'},
  {f: 'R', o: 'mod', a: [-81, 0, 0, 0, 1], b: [27, 0, 3], r: '0'}, // 0 is correct
  {f: 'Q', o: 'integrate', a: '3x^2', b: 0, r: '3x^2'},
  {f: 'Q', o: 'integrate', a: '3x^2', b: 1, r: 'x^3'},
  {f: 'Q', o: 'integrate', a: '3x^2', b: 2, r: '0.25x^4'},
  {f: 'R', o: 'add', a: "3x^2", b: "-x^2", r: "2x^2"},
  {f: 'Q', o: 'div', a: "3/2x^2-4x", b: "5x", r: "0.3x-0.8"},
  {f: 'Q', o: 'mod', a: "3/2x^2-4x", b: "5x", r: "0"},
  {f: 'Q', o: 'mod', a: "3x^2-4x", b: "5x", r: "0"},
  {f: 'Q', o: 'mod', a: "3x^2+2x+1", b: "x^2+1", r: "2x-2"},
  {f: 'Q', o: 'mod', a: "3x^2+2x+1", b: "2", r: "0"},
  {f: 'C', o: 'div', a: "25x^9", b: "5ix", r: "-5ix^8"},
  {f: 'C', o: 'div', a: "25x^9", b: "5x", r: "5x^8"},
  {f: 'C', o: 'mod', a: "25x^9", b: "5ix", r: "0"},
  {f: 'R', o: "mod", a: "5x^3+8x", b: "3x^2", r: "8x"},
  {f: 'Z7', o: 'pow', a: "9x^2+4", b: 3, r: "x^6+6x^4+5x^2+1"},
  {f: 'R', o: 'mul', a: "3x^3-1", b: 4, r: "12x^3-4"},
  {f: 'Q', o: 'integrate', a: '3x^2', b: 3, r: '0.05x^5'},
  {f: 'Q', o: 'integrate', a: '3x^2+3x^4', b: 2, r: '0.1x^6+0.25x^4'},
  {f: 'Q', o: 'derive', a: '3x^2', b: 0, r: '3x^2'},
  {f: 'Q', o: 'derive', a: '3x^2', b: 1, r: '6x'},
  {f: 'Q', o: 'derive', a: '3x^2', b: 2, r: '6'},
  {f: 'R', o: 'derive', a: '3x^2', b: 3, r: '0'},
  {f: 'R', o: 'derive', a: '3x^2', b: 4, r: '0'},
  {f: 'Q', o: 'derive', a: '3x^2', b: 5, r: '0'},
  {f: 'Q', o: 'pow', a: '3x', b: 9, r: '19683x^9'},
  {f: 'C', o: 'add', a: '3x+ix', b: "3ix^2+x^2", r: 'x^2+3ix^2+3x+ix'},
  {f: 'Q', o: 'div', a: '3x^4+2x^3+x^2-4', b: '2x^2-4', r: '1.5x^2+x+3.5'},
  {f: 'Q', o: 'mod', a: '3x^4+2x^3+x^2-4', b: '2x^2-4', r: '4x+10'},
  {f: 'Q', o: 'pow', a: '3x^2', b: 9, r: '19683x^18'},
  {f: 'Z7', o: 'pow', a: '3x^2', b: 9, r: '6x^18'},
  {f: 'C', o: 'pow', a: '3x^2+i', b: 2, r: '9x^4+6ix^2-1'},
  {f: 'R', o: 'add', a: 'x', b: "x^2", r: 'x^2+x'},
  {f: 'Z11', o: 'mod', a: '3x^2+x+7', b: '3x^2+x+7', r: '0'},
  {f: 'R', o: 'div', a: 'x^4+x^3+2x^2+1', b: 'x^3+2x^2+2', r: 'x-1'},
  {f: 'R', o: 'mod', a: 'x^4+x^3+2x^2+1', b: 'x^3+2x^2+2', r: '4x^2-2x+3'},
  {f: 'Q', o: "div", a: "x^4+x^3+2x^2+1", b: "2", r: "0.5x^4+0.5x^3+x^2+0.5"},
  {f: 'Q', o: "mod", a: "x^4+x^3+2x^2+1", b: "2", r: "0"},
  {f: 'R', o: 'add', a: '3x^3+4+6x^2+x', b: "-4", r: '3x^3+6x^2+x'},
  {f: 'R', o: 'add', a: '3x', b: "3x^3", r: '3x^3+3x'},
  {f: 'Z7', o: 'div', a: 'x^4+3x^3+6x+2', b: '2x^2+x+3', r: '4x^2+3x+3'},
  {f: 'Z7', o: 'mod', a: 'x^4+3x^3+6x+2', b: '2x^2+x+3', r: 'x'},
  {f: 'Q', o: 'derive', a: '-1/2x^2', b: 5, r: '0'},
  {f: 'Q', o: 'add', a: "-3x^6-13x^5-19x^4-23x^3-10x^2-10x+6", b: 0, r: '-3x^6-13x^5-19x^4-23x^3-10x^2-10x+6'},
  {f: 'Q', o: 'sub', a: "-18x^4-42x^3-48x^2-12x+12", b: 5, r: '-18x^4-42x^3-48x^2-12x+7'},
  {f: 'Q', o: 'mul', a: "55+3x^3+6x^5-x+4-1/2x^2+3x-x^2+x-4/8x^3", b: 1, r: '6x^5+2.5x^3-1.5x^2+3x+59'},
  {f: 'Q', o: 'mod', a: 'x^4+x^3+2x+2', b: 'x^2+2x+2', r: '4x+2'},
  {f: 'Q', o: 'div', a: 'x^4+x^3+2x+2', b: 'x^4+x^3+2x+2', r: '1'},
  {f: 'C', o: 'add', a: '1', b: '2', r: '3'},
  {f: 'C', o: 'add', a: '2i', b: 'i3', r: '5i'},
  {f: 'C', o: 'mul', a: '1-i', b: 'i+1', r: '2'},
  {f: 'C', o: 'add', a: 'x^2-ix^3', b: '12', r: '-ix^3+x^2+12'},
  {f: 'C', o: 'add', a: '-12', b: '12', r: '0'},
  {f: 'C', o: 'pow', a: '-i', b: 2, r: '-1'},
  {f: 'C', o: 'div', a: 'x-2', b: 'x-1', r: '1'},
  {f: 'C', o: 'div', a: 'x^2-2x+1', b: 'x', r: 'x-2'},
  {f: 'C', o: 'div', a: 'x^3-x^2+x-1', b: 'x-1', r: 'x^2+1'},
  {f: 'C', o: 'div', a: 'x^4-3x^3+3x^2-3x+2', b: 'x-1', r: 'x^3-2x^2+x-2'},
  {f: 'C', o: 'div', a: 'x^4-10x^3+26x^2-50x+105', b: 'x^2+5', r: 'x^2-10x+21'},
  {f: 'C', o: 'div', a: 'x^3-6x^2+11x-6', b: 'x-4', r: 'x^2-2x+3'},
  {f: 'C', o: 'div', a: 'x^5-5x^3+4x^2+x', b: 'x^2+3x+1', r: 'x^3-3x^2+3x-2'},
  {f: 'C', o: 'div', a: 'x^4-3x^2+3x+2', b: 'x^2-1', r: 'x^2-2'},
  {f: 'Q', o: 'gcd', a: '30x^2-90x+60', b: '50x^2-100x+50', r: 'x-1'},
  {f: 'Z5', o: 'gcd', a: '3x^3+4x^2+3', b: '3x^3+4x^2+3x+4', r: 'x+2'},
  {f: 'Z11', o: 'gcd', a: '3x^2+x+7', b: '3x^2+x+7', r: 'x^2+4x+6'},
  {f: 'Q', o: 'gcd', a: '3x^2+x+7', b: '3x^2+x+7', r: 'x^2+0.(3)x+2.(3)'},
  {f: 'Q', o: 'gcd', a: 'x^3+x^2+2x', b: '1+x^2', r: '1'},
  {f: 'Q', o: 'gcd', a: '15x^0', b: '3x^0', r: '1'},
  {f: 'Q', o: 'gcd', a: '30x^2-90x+60', b: '50x^2-100x+50', r: 'x-1'},
  {f: 'R', o: 'gcd', a: 'x^4+x^3+2x^2+1', b: 'x^3+2x^2+2', r: '1'},
  {f: 'C', o: 'gcd', a: "25x^9", b: "5ix", r: "x"},
  {f: 'Q', o: 'gcd', a: "x^4-9x^2-4x+12", b: "x^3+5x^2+2x-8", r: "x^2+x-2"},
  {f: 'Q', o: 'gcd', a: "30", b: "5", r: "1"},
  {f: 'Q', o: 'gcd', a: 'x^2+7x+6', b: 'x^2-5x-6', r: 'x+1'},
  {f: 'R', o: 'gcd', a: 'x^5+1', b: 'x^3+1', r: 'x+1'},
  {f: 'C', o: 'gcd', a: '-12', b: '0', r: '1'},
  {f: 'Q', o: 'gcd', a: '0', b: '0', r: '0'},
  {f: 'Q', o: 'gcd', a: 'x^2-4', b: 'x+2', r: 'x+2'},
  {f: 'Q', o: 'gcd', a: 'x^2+2', b: 'x^2-2', r: '1'},
  {f: 'Q', o: 'gcd', a: 'x^3-1', b: 'x^2-1', r: 'x-1'},
  {f: 'Q', o: 'gcd', a: 'x^4-5x^3+8x^2-10x+12', b: 'x^4+x^2-2', r: 'x^2+2'},
  {f: 'Q', o: 'gcd', a: 'x^4-16', b: 'x^4+8x^3+24x^2+32x+16', r: 'x+2'},
  {f: 'Q', o: 'gcd', a: 'x^5+x^4+x^3+x^2+x+1', b: 'x^3+x^2+x+1', r: 'x+1'},
  {f: 'Q', o: 'gcd', a: 'x^7+x^4+x^2+x^2+x+1', b: 'x^5+x^4+x^3+x^2+x+1', r: '1'},
  {f: 'Q', o: 'gcd', a: 'x^4+2x^3-2x^2+2x-3', b: 'x^4+x^3+2x^2+x+1', r: 'x^2+1'},
  {f: 'Q', o: 'gcd', a: 'x^3+x^2-3x-3', b: 'x^2+3x+2', r: 'x+1'},
  {f: 'Q', o: 'gcd', a: '12x^2+12x', b: '16x^2', r: 'x'},
  {f: 'C', o: 'gcd', a: 'x^2-1', b: 'x^2-3x+2', r: 'x-1'},
  {f: 'C', o: 'div', a: 'x^2-3', b: 'x+2', r: 'x-2'},
  {f: 'C', o: 'mod', a: 'x^2+3x+2', b: 'x+1', r: '0'},
  {f: 'C', o: 'mod', a: 'x^3+x^2-3x-3', b: 'x+1', r: '0'},
  {f: 'Z5', o: 'mod', a: '3x^3+4x^2+3', b: 'x+2', r: '0'},
  {f: 'Z5', o: 'mod', a: '3x^3+4x^2+3', b: '2x+4', r: '0'},
  {f: 'Z5', o: 'mod', a: '3x^3+4x^2+3x+4', b: 'x+2', r: '0'},
  {f: 'Z5', o: 'mod', a: '3x^3+4x^2+3x+4', b: '2x+4', r: '0'},
  {f: 'Z11', o: 'mod', a: '3x^2+x+7', b: 'x^2+4x+6', r: '0'},
  {f: 'Z11', o: 'mod', a: '3x^2+x+7', b: '3x^2+x+7', r: '0'},
  {f: 'Z11', o: 'mod', a: '3x^2+x+7', b: 'x^2+4x+6', r: '0'},
  {f: 'Z11', o: 'mod', a: '3x^2+x+7', b: '3x^2+x+7', r: '0'},
  {f: 'R', o: 'lc', a: '1', b: null, r: "1"},
  {f: 'R', o: 'lc', a: '4x+3', b: null, r: "4"},
  {f: 'R', o: 'lc', a: '5x^2+4x+3', b: null, r: "5"},
  {f: 'R', o: 'lm', a: '1', b: null, r: "1"},
  {f: 'R', o: 'lm', a: '4x+3', b: null, r: "4x"},
  {f: 'R', o: 'lm', a: '5x^2+4x+3', b: null, r: "5x^2"},
  {f: 'R', o: 'monic', a: 'x^2+10x+20', b: null, r: "x^2+10x+20"},
  {f: 'R', o: 'monic', a: '5x^2+10x+20', b: null, r: "x^2+2x+4"},
  {f: 'R', o: 'monic', a: '56', b: null, r: "1"},
  {f: 'Q', o: 'eval', a: '4/3x^2-3x+9', b: '2/3', r: '7.(592)'},
  {f: 'Q', o: 'eval', a: '4/3x^2-3x+9', b: '0', r: '9'},
  {f: 'R', o: 'eval', a: '0', b: 1, r: '0'},
  {f: 'R', o: 'eval', a: 'x-x', b: 199, r: '0'},
  {f: 'R', o: 'degree', a: 123, b: null, r: '0'},
  {f: 'R', o: 'degree', a: 'x^88+3', b: null, r: '88'},
  {f: 'R', o: 'add', a: '2', b: '2/1', r: '4'},
  {f: 'R', o: 'add', a: '2*2/2', b: '-1', r: '1'},
  {f: 'R', o: 'add', a: 'x', b: '3', r: 'x+3'},
  {f: 'R', o: 'add', a: '2x', b: '2x', r: '4x'},
  {f: 'R', o: 'add', a: '2*x', b: '2*x', r: '4x'},
  {f: 'R', o: 'add', a: '2*x^4', b: '2*x', r: '2x^4+2x'},
  {f: 'Q', o: 'add', a: 'x^8+1/2x', b: '1', r: 'x^8+0.5x+1'},
  {f: 'C', o: 'add', a: "2i*3*4x^8+1/2*x", b: '-0.5x', r: '24ix^8'},
  {f: 'Q', o: 'add', a: "3/4x^3+2/2x^2-1/3x+1/2", b: 0, r: "0.75x^3+x^2-0.(3)x+0.5"},
  {f: 'Q', o: 'add', a: ["1/2", "-1/3", "2/2", "3/4"], b: 0, r: "0.75x^3+x^2-0.(3)x+0.5"}
];



describe('Polynomial', function() {

  for (let i = 0; i < tests.length; i++) {

    (function(i) {

      const t = tests[i];

      it(t.a + " " + t.o + " " + t.b + " (in " + t.f + ")", function() {

        Poly.setField(t.f);
        let str;
        try {
          str = new Poly(t.a)[t.o](t.b, t.c).toString();
        } catch (e) {
          str = e;
        }
        assert.equal(t.r, str);
      });

    })(i);
  }
});

describe('Polynomial reciprocal', function() {

  it("Reciprocal polynomial of a complex polynomial with all coefficients", function() {
    Poly.setField("C");
      assert.equal(new Poly("2.x^4-5.ix^3+0.3x^2-0.1ix+4").reciprocal().toString(),
        new Poly("4x^4-0.1ix^3+0.3x^2-5.ix+2").toString());
  });

  it("Reciprocal polynomial of a real polynomial with vanishing lower terms", function() {
    Poly.setField("R");
      assert.equal(new Poly("2.x^4-5.x^3+0.3x^2").reciprocal().toString(),
        new Poly("0.3x^2-5.x+2").toString());
  });

  it('is an involution for polynomials with non-zero constant term', function () {
    Poly.setField('Z5');
    const polynomial = new Poly([1, 2, 3]);

    assert.equal(polynomial.reciprocal().toString(), 'x^2+2x+3');
    assert.equal(polynomial.reciprocal().reciprocal().equals(polynomial), true);
    assert.equal(new Poly(0).reciprocal().isZero(), true);
  });

});

describe('Polynomial fromRoots', function() {

  it("(x-1)(x-2i)(x-3) expands to x^3-4x^2-2ix^2+3x+8ix-6i", function() {
    Poly.setField("C");
    assert.equal(Poly.fromRoots([1., { 're':0. , 'im':2. }, 3.]).toString(),
                 "x^3-4x^2-2ix^2+3x+8ix-6i");
  });
  it("Depressed polynomial, (x-1)(x-2i)(x-3)x^4 expands to x^7-4x^6-2ix^6+3x^5+8ix^5-6ix^4", function() {
    Poly.setField("C");
    assert.equal(Poly.fromRoots([1., { 're':0. , 'im':2. }, 3., 0., 0., 0., 0.]).toString(),
                 "x^7-4x^6-2ix^6+3x^5+8ix^5-6ix^4");
  });

});

describe('Polynomial Horner', function() {

  it("1/2x^4-9x^2-4x+12 in horner", function() {

    Poly.setField("Q");
    assert.equal(new Poly("1/2x^4-9x^2-4x+12").toHorner(), "((0.5x^2-9)x-4)x+12");

    Poly.setField("C");
    assert.equal(new Poly("-2x^4-9x^2-4x+12i+3ix^2").toHorner(), "((-2x^2+(-9 + 3i))x-4)x+12i");
    assert.equal(new Poly("-ix^4+2x^4-9x^2-4x+12i+3ix^2").toHorner(), "(((2 - i)x^2+(-9 + 3i))x-4)x+12i");
    assert.equal(new Poly("x^4+2.ix^3+2x^2+6x").toHorner(), '(((1x+2i)x+2)x+6)x');
    assert.equal(new Poly("x^4+1.ix^3+2x^2+6x").toHorner(), '(((1x+i)x+2)x+6)x');

    Poly.setField("R");
    assert.equal(new Poly("3x^7-9x^4-4x^2+12").toHorner(), "((3x^3-9)x^2-4)x^2+12");
    assert.equal(new Poly("7x^2+x").toHorner(), "(7x+1)x");
    assert.equal(new Poly("7x^13+8x^7+3x^3+4x+2").toHorner(), "(((7x^6+8)x^4+3)x^2+4)x+2");

  });

  it("1/2x^4-9x^2-4x+12 in latex", function() {

    Poly.setField("Q");
    assert.equal(new Poly("1/2x^4-9x^2-4x+12").integrate().toLatex(), "\\frac{1}{10}x^5-3x^3-2x^2+12x");
  });
});

describe('Polynomial class and fields', function () {
  it('requires class construction with new', function () {
    assert.throws(() => Poly('x+1'), TypeError);
  });

  it('keeps the selected field bound to existing instances', function () {
    Poly.setField('Q');
    const rational = new Poly('1/3x');

    Poly.setField('R');

    assert.equal(rational.add('1/3x').toString(), '0.(6)x');
  });

  it('supports an isolated custom field', function () {
    const mod5 = createModuloField(5);
    const polynomial = new Poly('3x+4', mod5);

    assert.equal(polynomial.add('4x+3').toString(), '2x+2');
    assert.equal(polynomial.eval(2), 0);
  });

  it('normalizes modular coefficients and ignores leading zeros', function () {
    Poly.setField('Z3');
    assert.equal(new Poly([5, 0, -1]).toString(), '2x^2+2');

    const polynomial = new Poly([1, 0, 1, 0]);
    assert.equal(polynomial.degree(), 2);
    assert.equal(polynomial.lc(), 1);
  });
});

describe('Polynomial compose', function () {
  it('composes polynomials using Horner evaluation', function () {
    const polynomial = new Poly('2x^3-3x+5');

    assert.equal(polynomial.compose('x^2+1').toString(), '2x^6+6x^4+3x^2+4');
    assert.equal(polynomial.compose(2).toString(), String(polynomial.eval(2)));
    assert.equal(new Poly(7).compose('x^3-1').toString(), '7');
    assert.equal(new Poly(0).compose('x^3-1').toString(), '0');
  });

  it('composes polynomials in the bound coefficient field', function () {
    const mod5 = createModuloField(5);
    const polynomial = new Poly('2x^2+3x+4', mod5);

    assert.equal(polynomial.compose('2x+1').toString(), '3x^2+4x+4');

    Poly.setField('Z5');
    assert.equal(new Poly('4x^2+3x+2').compose('3x+1').toString(), 'x^2+3x+4');
  });
});

describe('Polynomial arithmetic edge cases', function () {
  it('handles zero and lower-degree dividends in polynomial division', function () {
    assert.equal(new Poly(0).div('x+1').toString(), '0');
    assert.equal(new Poly(0).mod('x+1').toString(), '0');
    assert.equal(new Poly(1).div('x^2+1').toString(), '0');
    assert.equal(new Poly(1).mod('x^2+1').toString(), '1');
    assert.throws(() => new Poly('x+1').div(0), RangeError);
    assert.throws(() => new Poly('x+1').mod(0), RangeError);
  });

  it('defines zero to the zeroth power as the multiplicative identity', function () {
    assert.equal(new Poly(0).pow(0).toString(), '1');
    assert.equal(new Poly(0).pow(5).toString(), '0');
    assert.equal(new Poly('x+1').pow(1).toString(), 'x+1');
    assert.throws(() => new Poly('x+1').pow(-1), RangeError);
  });

  it('applies formal differentiation in positive characteristic', function () {
    Poly.setField('Z5');

    assert.equal(new Poly('x^5').derive().toString(), '0');
    assert.equal(new Poly('2x^2+3x+4').derive().toString(), '4x+3');
  });
});

describe('Polynomial helpers', function () {
  it('compares normalized coefficients in the bound field', function () {
    assert.equal(new Poly('x^2+1').equals({ 0: 1, 2: 1 }), true);
    assert.equal(new Poly('x^2+1').equals('x^2-1'), false);

    Poly.setField('Z5');
    assert.equal(new Poly('x+1').equals('6x+6'), true);
  });

  it('identifies zero and monic polynomials', function () {
    const zero = new Poly('x-x');
    assert.equal(zero.isZero(), true);
    assert.equal(zero.degree(), -Infinity);
    assert.equal(zero.lc(), undefined);
    assert.equal(new Poly('x^2+2').isZero(), false);
    assert.equal(new Poly('x^2+2').isMonic(), true);
    assert.equal(new Poly('2x^2+2').isMonic(), false);
    assert.equal(new Poly(0).isMonic(), false);
  });

  it('tests polynomial divisibility', function () {
    assert.equal(new Poly('x^2-1').isDivisibleBy('x-1'), true);
    assert.equal(new Poly('x^2-1').isDivisibleBy('x+2'), false);
    assert.throws(() => new Poly('x^2-1').isDivisibleBy(0), RangeError);
  });

  it('returns the monic least common multiple', function () {
    assert.equal(
      new Poly('x^2-1').lcm('x^2+x-2').toString(),
      'x^3+2x^2-x-2',
    );
    assert.equal(new Poly('x^2-1').lcm(0).toString(), '0');

    Poly.setField('Z3');
    const left = new Poly([2, 2, 1, 1]);
    const right = new Poly([1, 2, 2, 1]);
    assert.equal(left.gcd(right).toString(), 'x^2+2');
    assert.equal(left.lcm(right).toString(), 'x^4+x^2+1');
  });

  it('exports dense coefficient arrays without truncation', function () {
    const polynomial = new Poly('3x^2+1');

    assert.deepEqual(polynomial.toArray(), [1, 0, 3]);
    assert.deepEqual(polynomial.toArray(5), [1, 0, 3, 0, 0]);
    assert.deepEqual(new Poly(0).toArray(), []);
    assert.deepEqual(new Poly(0).toArray(2), [0, 0]);
    assert.throws(() => polynomial.toArray(-1), RangeError);
    assert.throws(() => polynomial.toArray(1.5), RangeError);
  });
});

describe('Polynomial quotient rings', function () {
  it('reduces arithmetic modulo the defining polynomial', function () {
    Poly.setField('Z2');
    const ring = Poly.quotient('x^3+x+1');
    const alpha = ring.create('x');

    assert.equal(ring.modulus.toString(), 'x^3+x+1');
    assert.equal(ring.create('x^3+x+1').isZero(), true);
    assert.equal(alpha.pow(3).toString(), 'x+1');
    assert.equal(alpha.pow(7).equals(ring.one()), true);
    assert.equal(alpha.inverse().toString(), 'x^2+1');
    assert.equal(alpha.mul(alpha.inverse()).equals(1), true);
    assert.equal(alpha.div(alpha).equals(1), true);
  });

  it('evaluates polynomials at quotient-ring elements', function () {
    Poly.setField('Z2');
    const ring = Poly.quotient('x^3+x+1');
    const alpha = ring.create('x');

    assert.equal(ring.evaluate('x^3+x+1', alpha).toString(), '0');
    assert.equal(ring.isRoot('x^3+x+1', alpha), true);
    assert.equal(ring.evaluate('x^2+x+1', alpha).toString(), 'x^2+x+1');
    assert.equal(ring.isRoot('x^2+x+1', alpha), false);
    assert.equal(ring.evaluate(0, alpha).toString(), '0');
    assert.equal(ring.evaluate(1, alpha).toString(), '1');
    assert.equal(ring.isRoot(alpha.minimalPolynomial(), alpha), true);
  });

  it('evaluates in general quotient rings and rejects foreign elements', function () {
    Poly.setField('Z5');
    const ring = Poly.quotient('x^2-1');
    assert.equal(ring.isRoot('x^2-1', ring.create('x')), true);

    const foreign = Poly.quotient('x^2+1').create('x');
    assert.throws(() => ring.evaluate('x', foreign), RangeError);
    assert.throws(() => ring.isRoot('x', foreign), RangeError);
  });

  it('supports the ring operations on canonical representatives', function () {
    const ring = Poly.quotient('x^2+1');
    const value = ring.create('x+1');

    assert.equal(value.add('x-1').toString(), '2x');
    assert.equal(value.sub('x-1').toString(), '2');
    assert.equal(value.mul('x-1').toString(), '-2');
    assert.equal(value.neg().toString(), '-x-1');
    assert.equal(value.pow(0).equals(ring.one()), true);
    assert.equal(ring.zero().isZero(), true);
  });

  it('rejects non-units and incompatible quotient rings', function () {
    Poly.setField('Z5');
    const reducible = Poly.quotient('x^2-1');
    assert.throws(() => reducible.create('x-1').inverse(), RangeError);
    assert.throws(() => reducible.zero().inverse(), RangeError);

    const left = Poly.quotient('x^2+1').create('x');
    const right = Poly.quotient('x^2+2').create('x');
    assert.throws(() => left.add(right), RangeError);
  });

  it('requires a non-constant modulus and preserves an instance field', function () {
    assert.throws(() => Poly.quotient(0), RangeError);
    assert.throws(() => Poly.quotient(1), RangeError);

    Poly.setField('Q');
    const modulus = new Poly('x^2+1');
    Poly.setField('R');
    const ring = Poly.quotient(modulus);
    assert.equal(ring.create('1/3x').add('1/3x').toString(), '0.(6)x');
  });

  it('computes orders and primitive elements in finite extension fields', function () {
    Poly.setField('Z2');
    const ring = Poly.quotient('x^3+x+1');
    const alpha = ring.create('x');

    assert.equal(ring.order(), 8);
    assert.equal(ring.isField(), true);
    assert.equal(ring.elements().length, 8);
    assert.equal(ring.primitiveElements().length, 6);
    assert.equal(alpha.additiveOrder(), 2);
    assert.equal(ring.zero().additiveOrder(), 1);
    assert.equal(alpha.multiplicativeOrder(), 7);
    assert.equal(alpha.isPrimitive(), true);
    assert.equal(ring.one().isPrimitive(), false);
    assert.throws(() => ring.zero().multiplicativeOrder(), RangeError);
  });

  it('enumerates generated additive and multiplicative subgroups', function () {
    Poly.setField('Z2');
    const ring = Poly.quotient('x^3+x+1');
    const alpha = ring.create('x');

    assert.deepEqual(alpha.additiveSubgroup().map(String), ['0', 'x']);
    assert.deepEqual(ring.zero().additiveSubgroup().map(String), ['0']);
    assert.deepEqual(alpha.multiplicativeSubgroup().map(String), [
      '1',
      'x',
      'x^2',
      'x+1',
      'x^2+x',
      'x^2+x+1',
      'x^2+1',
    ]);
    assert.deepEqual(ring.one().multiplicativeSubgroup().map(String), ['1']);
    assert.throws(() => alpha.additiveSubgroup(1), RangeError);
    assert.throws(() => alpha.multiplicativeSubgroup(6), RangeError);
    assert.throws(() => alpha.additiveSubgroup(-1), RangeError);
    assert.throws(() => ring.zero().multiplicativeSubgroup(), RangeError);
  });

  it('derives minimal polynomials from Frobenius orbits', function () {
    Poly.setField('Z2');
    const ring = Poly.quotient('x^4+x+1');
    const alpha = ring.create('x');
    const subfieldElement = alpha.pow(5);

    assert.equal(alpha.minimalPolynomial().toString(), 'x^4+x+1');
    assert.equal(subfieldElement.multiplicativeOrder(), 3);
    assert.equal(subfieldElement.minimalPolynomial().toString(), 'x^2+x+1');
    assert.equal(ring.zero().minimalPolynomial().toString(), 'x');
  });

  it('enumerates distinct primitive polynomials', function () {
    Poly.setField('Z2');
    const cubicField = Poly.quotient('x^3+x+1');
    const quarticField = Poly.quotient('x^4+x+1');

    assert.deepEqual(cubicField.primitivePolynomials().map(String), [
      'x^3+x+1',
      'x^3+x^2+1',
    ]);
    assert.deepEqual(quarticField.primitivePolynomials().map(String), [
      'x^4+x+1',
      'x^4+x^3+1',
    ]);
    assert.equal(
      cubicField.primitiveElements().length,
      cubicField.primitivePolynomials().length * cubicField.modulus.degree(),
    );
    for (const polynomial of quarticField.primitivePolynomials()) {
      assert.equal(polynomial.isIrreducible(), true);
      assert.equal(polynomial.degree(), 4);
    }
  });

  it('guards finite-field operations and enumeration size', function () {
    Poly.setField('Z5');
    const reducible = Poly.quotient('x^2-1');
    assert.equal(reducible.isField(), false);
    assert.throws(() => reducible.primitiveElements(), RangeError);
    assert.throws(() => reducible.primitivePolynomials(), RangeError);
    assert.throws(() => reducible.create('x').minimalPolynomial(), RangeError);
    assert.throws(() => Poly.quotient('x^4+2').elements(100), RangeError);
    assert.throws(() => Poly.quotient('x^4+2').primitivePolynomials(100), RangeError);

    Poly.setField('R');
    assert.throws(() => Poly.quotient('x^2+1').order(), TypeError);
  });
});

describe('Polynomial cyclic codes', function () {
  it('builds a shifted basis from a generator polynomial', function () {
    Poly.setField('Z2');
    const generator = new Poly('x^4+x^3+x^2+1');

    assert.deepEqual(generator.cyclicCodeBasis(7).map(String), [
      'x^4+x^3+x^2+1',
      'x^5+x^4+x^3+x',
      'x^6+x^5+x^4+x^2',
    ]);
  });

  it('normalizes generators and validates the code contract', function () {
    Poly.setField('Z3');
    assert.deepEqual(new Poly('2x+2').cyclicCodeBasis(2).map(String), ['x+1']);
    assert.throws(() => new Poly('x+1').cyclicCodeBasis(0), RangeError);
    assert.throws(() => new Poly(0).cyclicCodeBasis(7), RangeError);
    assert.throws(() => new Poly('x^2+x+1').cyclicCodeBasis(4), RangeError);

    Poly.setField('R');
    assert.throws(() => new Poly('x+1').cyclicCodeBasis(2), TypeError);
  });

  it('finds every generator polynomial for a code length and dimension', function () {
    Poly.setField('Z2');

    assert.deepEqual(Poly.cyclicCodeGenerators(7, 4).map(String), [
      'x^3+x^2+1',
      'x^3+x+1',
    ]);
    assert.deepEqual(Poly.cyclicCodeGenerators(7, 3).map(String), [
      'x^4+x^2+x+1',
      'x^4+x^3+x^2+1',
    ]);
    assert.deepEqual(Poly.cyclicCodeGenerators(4, 2).map(String), ['x^2+1']);
    assert.deepEqual(Poly.cyclicCodeGenerators(7, 5), []);
  });

  it('supports zero and full dimensions and bounds generator searches', function () {
    Poly.setField('Z2');
    assert.deepEqual(Poly.cyclicCodeGenerators(3, 0).map(String), ['x^3+1']);
    assert.deepEqual(Poly.cyclicCodeGenerators(3, 3).map(String), ['1']);
    assert.throws(() => Poly.cyclicCodeGenerators(7, 4, 1), RangeError);
    assert.throws(() => Poly.cyclicCodeGenerators(0, 0), RangeError);
    assert.throws(() => Poly.cyclicCodeGenerators(7, -1), RangeError);
    assert.throws(() => Poly.cyclicCodeGenerators(7, 8), RangeError);

    Poly.setField('R');
    assert.throws(() => Poly.cyclicCodeGenerators(7, 4), TypeError);
  });
});

describe('Polynomial finite-field enumeration', function () {
  it('enumerates polynomials in base-field coefficient order', function () {
    Poly.setField('Z2');
    assert.deepEqual(Poly.enumerate(2).map(String), [
      '0',
      '1',
      'x',
      'x+1',
      'x^2',
      'x^2+1',
      'x^2+x',
      'x^2+x+1',
    ]);
  });

  it('compares polynomials and carries into the next degree', function () {
    Poly.setField('Z3');
    assert.equal(new Poly('2x+2').successor().toString(), 'x^2');
    assert.equal(new Poly('2x+1').compare('2x+2'), -1);
    assert.equal(new Poly('2x+2').compare('2x+2'), 0);
    assert.equal(new Poly('x^2').compare('2x+2'), 1);
  });

  it('validates finite enumeration bounds and fields', function () {
    Poly.setField('Z2');
    assert.throws(() => Poly.enumerate(-1), RangeError);
    assert.throws(() => Poly.enumerate(10, 100), RangeError);

    Poly.setField('R');
    assert.throws(() => Poly.enumerate(2), TypeError);
    assert.throws(() => new Poly('x').successor(), TypeError);
  });
});

describe('Polynomial finite-field factorization', function () {
  it('factors square-free polynomials over prime fields', function () {
    Poly.setField('Z2');
    const polynomial = new Poly('x^5+x^2+x+1');
    const factors = polynomial.factor();

    assert.deepEqual(factors.map(String), ['x+1', 'x+1', 'x^3+x+1']);
    assert.equal(
      factors.reduce((product, factor) => product.mul(factor), new Poly(1)).equals(polynomial),
      true,
    );
  });

  it('preserves multiplicities in inseparable polynomials', function () {
    Poly.setField('Z2');
    assert.deepEqual(
      new Poly('x^4+x^2+1').factor().map(String),
      ['x^2+x+1', 'x^2+x+1'],
    );

    Poly.setField('Z3');
    assert.deepEqual(
      new Poly('x^6+2x^3+1').factor().map(String),
      Array(6).fill('x+1'),
    );
  });

  it('identifies irreducible and reducible polynomials', function () {
    Poly.setField('Z3');

    assert.equal(new Poly('x^5+2x+2').isIrreducible(), true);
    assert.equal(new Poly('x^5-1').isIrreducible(), false);
    assert.equal(new Poly('2x+2').isIrreducible(), true);
    assert.equal(new Poly(2).isIrreducible(), false);
  });

  it('retains the leading unit of non-monic inputs', function () {
    Poly.setField('Z3');

    assert.deepEqual(new Poly('2x+2').factor().map(String), ['2', 'x+1']);
    assert.deepEqual(new Poly(2).factor().map(String), ['2']);
    assert.deepEqual(new Poly(1).factor(), []);
  });

  it('rejects zero and unsupported coefficient fields', function () {
    Poly.setField('Z2');
    assert.throws(() => new Poly(0).factor(), RangeError);

    Poly.setField('R');
    assert.throws(() => new Poly('x^2-1').factor(), TypeError);
    assert.throws(() => new Poly('x^2-1').isIrreducible(), TypeError);

    Poly.setField('Z4');
    assert.throws(() => new Poly('x^2+1').factor(), TypeError);
  });
});

describe('Polynomial dispersion', function () {
  it('finds the greatest integral distance between polynomial roots', function () {
    Poly.setField('Q');

    assert.equal(new Poly('x^3-7x^2+10x').dispersion(), 5);
    assert.equal(new Poly('x^2-25x').dispersion(), 25);
    assert.equal(new Poly('x^2+1').dispersion(), 0);
    assert.equal(new Poly(7).dispersion(), 0);

    Poly.setField('C');
    assert.equal(new Poly('x^2-3x').dispersion(), 3);

    Poly.setField('R');
    assert.equal(new Poly('x^2-2.2x+0.21').dispersion(), 2);
  });

  it('handles bounded custom and positive-characteristic fields', function () {
    const mod5 = createModuloField(5);
    const polynomial = new Poly('x^2-2x', mod5);

    assert.throws(() => polynomial.dispersion(), TypeError);
    assert.equal(polynomial.dispersion(2), 2);

    Poly.setField('Z5');
    assert.equal(new Poly('x^2-2x').dispersion(), Infinity);
  });

  it('validates an explicit dispersion bound', function () {
    const polynomial = new Poly('x^2-x');

    assert.throws(() => polynomial.dispersion(-1), RangeError);
    assert.throws(() => polynomial.dispersion(1.5), RangeError);
  });
});

describe('Polynomial square root', function () {
  it('extracts the polynomial square roots from the documented examples', function () {
    const first = new Poly('9x^4+6x^3-11x^2-4x+4');
    const firstRoot = first.sqrt();
    assert.ok(firstRoot instanceof Poly);
    assert.equal(firstRoot.toString(), '3x^2+x-2');
    assert.equal(firstRoot.mul(firstRoot).toString(), first.toString());

    const second = new Poly('x^6-2x^5-x^4+3x^2+2x+1');
    const secondRoot = second.sqrt();
    assert.ok(secondRoot instanceof Poly);
    assert.equal(secondRoot.toString(), 'x^3-x^2-x-1');
    assert.equal(secondRoot.mul(secondRoot).toString(), second.toString());
  });

  it('extracts exact rational and complex square roots', function () {
    Poly.setField('Q');
    const rationalRoot = new Poly('3/2x^2-2/3x+5/7');
    const rationalSquare = rationalRoot.mul(rationalRoot);
    assert.equal(rationalSquare.sqrt()?.toString(), rationalRoot.toString());

    Poly.setField('C');
    const complexSquare = new Poly('-x^2');
    const complexRoot = complexSquare.sqrt();
    assert.equal(complexRoot?.toString(), 'ix');
    assert.equal(complexRoot?.mul(complexRoot).toString(), complexSquare.toString());
  });

  it('supports irrational and decimal roots over the real field', function () {
    const irrational = new Poly('2x^2').sqrt();
    assert.ok(irrational instanceof Poly);
    assert.ok(Math.abs(irrational.coeff[1] - Math.SQRT2) < Number.EPSILON);

    const decimalSquare = new Poly('0.01x^2+0.04x+0.04');
    assert.equal(decimalSquare.sqrt()?.toString(), '0.1x+0.2');
  });

  it('extracts square roots over prime fields', function () {
    Poly.setField('Z7');
    const oddCharacteristic = new Poly('2x^2+5x+4');
    const oddRoot = oddCharacteristic.sqrt();
    assert.equal(oddRoot?.toString(), '3x+2');
    assert.equal(oddRoot?.mul(oddRoot).toString(), oddCharacteristic.toString());

    Poly.setField('Z5');
    const tonelliShanks = new Poly('4x^2+4x+1');
    assert.equal(tonelliShanks.sqrt()?.toString(), '2x+1');

    Poly.setField('Z2');
    const characteristicTwo = new Poly('x^4+x^2+1');
    const binaryRoot = characteristicTwo.sqrt();
    assert.equal(binaryRoot?.toString(), 'x^2+x+1');
    assert.equal(binaryRoot?.mul(binaryRoot).toString(), characteristicTwo.toString());
  });

  it('returns null when no polynomial square root exists', function () {
    Poly.setField('Q');
    assert.equal(new Poly('x^2+1').sqrt(), null);
    assert.equal(new Poly('x^3+1').sqrt(), null);
    assert.equal(new Poly('2x^2').sqrt(), null);

    Poly.setField('R');
    assert.equal(new Poly('-x^2').sqrt(), null);
  });

  it('handles zero and reports unsupported coefficient fields', function () {
    assert.equal(new Poly(0).sqrt()?.toString(), '0');

    const mod5 = createModuloField(5);
    assert.throws(
      () => new Poly('x^2+2x+1', mod5).sqrt(),
      TypeError,
    );

    Poly.setField('H');
    assert.throws(() => new Poly('x^2').sqrt(), TypeError);
  });

  it('uses a custom coefficient square root when provided', function () {
    const field = {
      ...createModuloField(5),
      characteristic: 5,
      sqrt: value => {
        for (let candidate = 0; candidate < 5; candidate += 1) {
          if (candidate * candidate % 5 === value) {
            return candidate;
          }
        }
        return null;
      },
    };
    const square = new Poly('4x^2+4x+1', field);

    assert.equal(square.sqrt()?.mul(square.sqrt()).toString(), square.toString());
  });
});

describe('Polynomial validation and tracing', function () {
  it('rejects malformed expressions and exponents', function () {
    assert.throws(() => new Poly('x^-1'), SyntaxError);
    assert.throws(() => new Poly('x^1.5'), SyntaxError);
    assert.throws(() => new Poly('x^^2'), SyntaxError);
    assert.throws(() => new Poly({ '-1': 2 }), RangeError);
    assert.throws(() => new Poly('x').pow(1.5), RangeError);
  });

  it('retains the result alias and division trace', function () {
    const polynomial = new Poly('x^2-1');
    assert.equal(polynomial.result(3), polynomial.eval(3));

    Poly.trace = true;
    assert.equal(polynomial.div('x-1').toString(), 'x+1');
    assert.ok(Array.isArray(Poly.trace));
    assert.ok(Poly.trace.every(step => step instanceof Poly));
    assert.equal(Poly.trace.at(-1).toString(), '0');
  });

  it('validates built-in and custom fields', function () {
    assert.throws(() => Poly.setField('Z1'), RangeError);
    assert.throws(() => Poly.setField({ add() {} }), TypeError);
  });
});

