# Polynomial.js

[![NPM Package](https://img.shields.io/npm/v/polynomial.svg?style=flat)](https://www.npmjs.com/package/polynomial)
[![MIT license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)

Polynomial.js is a TypeScript library for parsing and computing univariate polynomials. It supports real, rational, complex, quaternion, modular, and custom coefficient fields.

The library provides immutable arithmetic operations, derivatives, integrals, Euclidean division, greatest common divisors, reciprocal polynomials, and Horner evaluation. It does not calculate roots, but it can construct a polynomial from known roots.

## Installation

```bash
npm install polynomial
```

Polynomial.js requires Node.js 20 or newer when used in Node.js.

## Quick Start

```js
import Polynomial from 'polynomial';

Polynomial.setField('R');

const polynomial = Polynomial.fromRoots([1, 2, 3]);

console.log(polynomial.toString());          // x^3-6x^2+11x-6
console.log(polynomial.eval(2));             // 0
console.log(polynomial.derive().toString()); // 3x^2-12x+11
```

Every polynomial is an instance of the `Polynomial` class and must be constructed with `new`. Arithmetic methods return new instances and do not mutate the receiver.

```js
const left = new Polynomial('3x^2');
const sum = left.add('-x^2');

console.log(left.toString()); // 3x^2
console.log(sum.toString());  // 2x^2
```

CommonJS is supported as well:

```js
const Polynomial = require('polynomial');
```

## Inputs

The constructor and arithmetic methods accept polynomial strings, numbers, coefficient arrays, coefficient maps, and other `Polynomial` instances.

```js
new Polynomial('23x^4+98x^2+4');
new Polynomial(55);
new Polynomial([1, 2, 3]);       // 3x^2+2x+1
new Polynomial({ 3: 4, 5: 9 });  // 9x^5+4x^3
```

String exponents must be non-negative integers. Whitespace is ignored. Coefficients may contain multiplication and division, while nested polynomial expressions and parentheses are not parsed.

```js
new Polynomial('2 * 3x^2 + 4/5x - 1');
```

The selected field parses each coefficient, enabling field-specific forms such as fractions, repeating decimals, and complex numbers.

```js
Polynomial.setField('Q');
console.log(new Polynomial('5/3x^3+4/3x').toString());

Polynomial.setField('C');
console.log(new Polynomial('23ix^4+98x^2+i').toString());
```

## Coefficient Fields

The built-in fields are:

| Name | Coefficients | Implementation |
| --- | --- | --- |
| `R` | Real numbers | JavaScript numbers |
| `Q` | Rational numbers | [Fraction.js](https://github.com/rawify/Fraction.js) |
| `C` | Complex numbers | [Complex.js](https://github.com/rawify/Complex.js) |
| `H` | Quaternions | [Quaternion.js](https://github.com/rawify/Quaternion.js) |
| `Zp` | Integers modulo `p` | Built in, for example `Z7` |

`Polynomial.setField()` selects the default field for future instances. Each instance captures that field when it is constructed, so changing the default does not alter existing polynomials or their subsequent operations.

```js
Polynomial.setField('Q');
const rational = new Polynomial('1/3x');

Polynomial.setField('R');
console.log(rational.add('1/3x').toString()); // 0.(6)x
```

Division in `Zp` requires the divisor to have a multiplicative inverse. Use a prime modulus when all non-zero coefficients must be invertible.

### Custom Fields

A custom field can be installed as the default with `setField()` or passed directly to the constructor. Passing it to the constructor avoids shared configuration and keeps the field local to that polynomial and all results derived from it.

```ts
import Polynomial, { type PolynomialField } from 'polynomial';

const field: PolynomialField<number> = {
  add: (left, right) => left + Number(right),
  sub: (left, right) => left - Number(right),
  mul: (left, right) => left * Number(right),
  div: (left, right) => left / Number(right),
  parse: value => Number(value),
  empty: value => value === undefined || value === 0,
  pow: (base, exponent) => base ** exponent,
  equals: (left, right) => left === Number(right),
  characteristic: 0,
  magnitude: value => Math.abs(value),
  sqrt: value => value < 0 ? null : Math.sqrt(value),
};

const polynomial = new Polynomial('x^2+2x+1', field);
```

## API

### Construction and Configuration

- `new Polynomial(input?, field?)` creates a polynomial.
- `Polynomial.fromRoots(roots, field?)` creates a monic polynomial with the given roots.
- `Polynomial.quotient(modulus, field?)` creates the quotient ring defined by a non-constant modulus polynomial.
- `Polynomial.enumerate(maxDegree, maxElements?)` or `Polynomial.enumerate(maxDegree, field, maxElements?)` enumerates polynomials over a finite prime field.
- `Polynomial.cyclicCodeGenerators(length, dimension, maxGenerators?)` or `Polynomial.cyclicCodeGenerators(length, dimension, field, maxGenerators?)` finds all generator polynomials for a cyclic code.
- `Polynomial.setField(field)` selects `R`, `Q`, `C`, `H`, `Zp`, or a custom field as the default for future instances.

### Arithmetic

- `add(value)` returns the sum.
- `sub(value)` returns the difference.
- `mul(value)` returns the product.
- `addmul(left, right)` adds the product of two values.
- `div(value)` returns the quotient of polynomial long division.
- `mod(value)` returns the remainder of polynomial long division.
- `gcd(value)` returns the monic greatest common divisor.
- `lcm(value)` returns the monic least common multiple.
- `isDivisibleBy(value)` reports whether polynomial division has zero remainder.
- `factor()` returns the leading unit followed by monic irreducible factors with multiplicity.
- `isIrreducible()` reports whether a positive-degree polynomial is irreducible.
- `cyclicCodeBasis(length)` returns the shifted basis generated by a divisor of $x^n-1$.
- `compare(value)` compares finite-field polynomials in coefficient order.
- `successor()` returns the next polynomial in that order.
- `pow(exponent)` raises the polynomial to a non-negative integer power.
- `neg()` negates every coefficient.
- `reciprocal()` reverses the coefficient order relative to the degree.
- `monic()` divides every coefficient by the leading coefficient.
- `dispersion(maxShift?)` returns the greatest non-negative integer `k` for which `gcd(p(x), p(x+k))` is non-constant.
- `sqrt()` returns a polynomial `g` satisfying `g² = p`, or `null` when no polynomial square root exists.

For `R`, `Q`, and `C`, `dispersion()` derives a finite search bound from the coefficients. A custom field can provide `magnitude()` and an optional `characteristic` to enable the same behavior. Otherwise, pass an explicit non-negative `maxShift`. Non-constant polynomials over a field with positive characteristic have infinite dispersion because shifts by multiples of the characteristic leave the polynomial unchanged.

```js
Polynomial.setField('Q');

const polynomial = new Polynomial('x^3-7x^2+10x'); // roots: 0, 2, 5
console.log(polynomial.dispersion());               // 5
console.log(polynomial.dispersion(4));              // 3
```

`sqrt()` uses the root with the principal leading coefficient. Its negation is the second square root. Exact roots are supported over `Q` and prime fields, while `R` also permits irrational coefficients and `C` uses the principal complex coefficient root. Custom fields must provide `sqrt(value)`. Quaternion coefficients are not supported because the coefficient recurrence assumes commutative multiplication.

```js
const square = new Polynomial('9x^4+6x^3-11x^2-4x+4');
const root = square.sqrt();

console.log(root?.toString());           // 3x^2+x-2
console.log(root?.mul(root).toString()); // 9x^4+6x^3-11x^2-4x+4
```

### Calculus and Evaluation

- `derive(order = 1)` returns the requested derivative.
- `integrate(order = 1)` returns the requested antiderivative with zero integration constants.
- `eval(value)` evaluates the polynomial using [Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method).
- `compose(value)` substitutes a polynomial for the indeterminate using Horner's method.
- `result(value)` is a deprecated alias for `eval(value)`.

```js
const polynomial = new Polynomial('2x^3-3x+5');
console.log(polynomial.compose('x^2+1').toString()); // 2x^6+6x^4+3x^2+4
```

### Quotient Rings

`Polynomial.quotient(modulus, field?)` constructs $K[x]/(m)$. Every element is reduced to its canonical representative of degree less than `modulus.degree()`. A ring created from an existing polynomial retains that polynomial's coefficient field.

```js
Polynomial.setField('Z2');

const ring = Polynomial.quotient('x^3+x+1');
const alpha = ring.create('x');

console.log(alpha.pow(3).toString());    // x+1
console.log(alpha.pow(7).toString());    // 1
console.log(alpha.inverse().toString()); // x^2+1
```

The ring provides `create(value)`, `zero()`, `one()`, and its monic `modulus`. `evaluate(polynomial, value)` evaluates a polynomial at a residue class using Horner's method, while `isRoot(polynomial, value)` tests whether that value is a root. Residue classes provide `add`, `sub`, `mul`, `div`, `neg`, `pow`, `inverse`, `equals`, and `isZero`. Division and inversion throw `RangeError` when the divisor is not a unit. Residues from different quotient-ring instances cannot be combined.

```js
Polynomial.setField('Z2');

const ring = Polynomial.quotient('x^3+x+1');
const alpha = ring.create('x');

console.log(ring.evaluate('x^2+x+1', alpha).toString()); // x^2+x+1
console.log(ring.isRoot('x^3+x+1', alpha));              // true
```

Unlike `Polynomial.eval()`, which returns a coefficient-field value, `ring.evaluate()` returns a residue class. Evaluation is defined for every quotient ring and does not require the modulus to be irreducible. Values already belonging to another quotient-ring instance are rejected.

When the coefficient field is `Zp` and the modulus is irreducible, the quotient is the finite field $\text{GF}(p^n)$. `isField()` checks this condition, and `order()` returns $p^n$. The ring can enumerate its `elements()`, `primitiveElements()`, and `primitivePolynomials()`. A residue class provides `additiveOrder()`, `additiveSubgroup()`, `multiplicativeOrder()`, `multiplicativeSubgroup()`, `isPrimitive()`, and `minimalPolynomial()`.

```js
Polynomial.setField('Z2');

const field = Polynomial.quotient('x^4+x+1');
const alpha = field.create('x');
const beta = alpha.pow(5);

console.log(field.order());                        // 16
console.log(alpha.multiplicativeOrder());          // 15
console.log(alpha.isPrimitive());                  // true
console.log(beta.multiplicativeSubgroup().map(String));
// ['1', 'x^2+x', 'x^2+x+1']
console.log(beta.minimalPolynomial().toString()); // x^2+x+1
console.log(field.primitivePolynomials().map(String));
// ['x^4+x+1', 'x^4+x^3+1']
```

`primitivePolynomials()` returns the distinct minimal polynomials of the primitive elements, ordered by finite-field coefficient order. Frobenius conjugates are deduplicated, so each primitive polynomial occurs exactly once.

The subgroup methods list the identity first, followed by successive multiples or powers of the residue. Their optional limit defaults to one million elements. `elements(maxElements?)`, `primitiveElements(maxElements?)`, and `primitivePolynomials(maxElements?)` use the same default and throw `RangeError` before enumeration if the field is larger. Multiplicative operations and minimal polynomials require an irreducible modulus; ordinary quotient-ring arithmetic does not.

### Finite-Field Factorization

`factor()` and `isIrreducible()` operate over prime coefficient fields `Zp`. Factorization uses square-free decomposition, Frobenius roots for inseparable polynomials, and Berlekamp's algorithm. Repeated factors occur repeatedly in the returned array.

```js
Polynomial.setField('Z2');

const polynomial = new Polynomial('x^5+x^2+x+1');
console.log(polynomial.factor().map(String));
// ['x+1', 'x+1', 'x^3+x+1']
```

For a non-monic polynomial, the leading unit is the first factor, so multiplying all returned factors reconstructs the input. The factorization of `1` is empty, while factoring the zero polynomial throws `RangeError`. Fields of characteristic zero and composite modular rings are rejected.

### Cyclic Codes

If $g$ divides $x^n-1$ over a prime field, `g.cyclicCodeBasis(n)` returns the standard basis $g,xg,\ldots,x^{k-1}g$, where $k=n-\deg g$.

```js
Polynomial.setField('Z2');

const generator = new Polynomial('x^4+x^3+x^2+1');
console.log(generator.cyclicCodeBasis(7).map(String));
// ['x^4+x^3+x^2+1', 'x^5+x^4+x^3+x', 'x^6+x^5+x^4+x^2']
```

The generator is normalized to monic form. Invalid lengths, the zero polynomial, and generators that do not divide $x^n-1$ are rejected.

`Polynomial.cyclicCodeGenerators(n, k)` returns every monic divisor of $x^n-1$ having degree $n-k$. Each result therefore generates a cyclic $[n,k]$ code and can be passed to `cyclicCodeBasis(n)`.

```js
Polynomial.setField('Z2');

const generators = Polynomial.cyclicCodeGenerators(7, 4);
console.log(generators.map(String));
// ['x^3+x^2+1', 'x^3+x+1']

console.log(generators[0].cyclicCodeBasis(7).length); // 4
```

Repeated irreducible factors are retained when the characteristic divides the code length. Dimensions from zero through $n$ are accepted: $k=0$ returns $x^n-1$, while $k=n$ returns $1$. An empty array means that no divisor has the required degree. The optional result limit defaults to one million generators.

### Finite-Field Enumeration

`compare()` and `successor()` use base-$p$ coefficient order with the constant coefficient as the least significant digit. `Polynomial.enumerate(d)` returns all $p^{d+1}$ polynomials of degree at most $d$ in that order.

```js
Polynomial.setField('Z2');

console.log(Polynomial.enumerate(1).map(String));
// ['0', '1', 'x', 'x+1']

console.log(new Polynomial('x+1').successor().toString()); // x^2
```

Enumeration defaults to at most one million results. Pass the limit as the second argument when using the active field, or as the third argument when providing a custom field.

### Inspection and Formatting

- `degree()` returns the polynomial degree, or `-Infinity` for the zero polynomial.
- `lc()` returns the leading coefficient, or `undefined` for the zero polynomial.
- `lm()` returns the leading monomial.
- `equals(value)` compares normalized coefficients in the current field.
- `isZero()` reports whether the polynomial is zero.
- `isMonic()` reports whether the leading coefficient is one.
- `toArray(width?)` returns coefficients in ascending exponent order, optionally padded with zeros.
- `clone()` returns an independent polynomial with the same coefficients and field.
- `toString()` returns the plain-text representation.
- `toLatex()` returns a LaTeX representation.
- `toHorner()` returns a Horner-scheme representation.
- `coeff` exposes the coefficient map indexed by exponent.

```js
const polynomial = new Polynomial('23x^4+98x^2+4');

console.log(polynomial.coeff[4]); // 23
console.log(polynomial.degree()); // 4
console.log(polynomial.toArray()); // [4, 0, 98, 0, 23]
```

### Division Trace

Set `Polynomial.trace` to `true` before division to retain the subtracted terms and final remainder. The property is replaced with an array of polynomial instances.

```js
Polynomial.trace = true;

new Polynomial('x^4+3x^3+2x^2+6x').div('x+3');
console.log(Polynomial.trace.map(step => step.toString()));
// ['x^4+3x^3', '2x^2+6x', '0']

Polynomial.trace = null;
```

## Browser Usage

The browser builds include all built-in coefficient implementations. No separate Fraction.js, Complex.js, or Quaternion.js scripts are required.

Use the global IIFE build directly:

```html
<script src="https://unpkg.com/polynomial@2.0.0/dist/polynomial.min.js"></script>
<script>
  Polynomial.setField('C');
  console.log(new Polynomial('x-i').mul('x+i').toString());
</script>
```

Or import the standalone browser ES module:

```js
import Polynomial from 'https://unpkg.com/polynomial@2.0.0/dist/polynomial.min.mjs';
```

AMD loaders receive the class as an anonymous module from `polynomial.min.js`.

## TypeScript

The package ships declarations for ESM and CommonJS. ESM exposes the class as both the default and named `Polynomial` export, together with `BuiltInField`, `PolynomialField`, `CoefficientMap`, `PolynomialInput`, `PolynomialQuotientRing`, `PolynomialResidue`, and `PolynomialResidueInput` types.

CommonJS namespace aliases are available as `Polynomial.FieldName`, `Polynomial.Field<C>`, `Polynomial.Coefficients<C>`, `Polynomial.Input<C>`, `Polynomial.QuotientRing<C>`, `Polynomial.Residue<C>`, and `Polynomial.ResidueInput<C>`.

## Errors

Malformed polynomial syntax throws `SyntaxError`. Invalid fields and exponents, division by zero, and non-invertible modular coefficients throw `RangeError` or `TypeError` as appropriate.

## Development

```bash
npm install
npm test
```

`npm test` performs a strict TypeScript build, checks ESM and CommonJS declaration consumers, and runs the Node.js and standalone browser tests.

## License

Copyright (c) 2026 [Robert Eisele](https://raw.org/)

Licensed under the MIT license.
