# Polynomial.js

[![NPM Package](https://img.shields.io/npm/v/polynomial.svg?style=flat)](https://npmjs.org/package/polynomial "View this project on npm")
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

Polynomials are defined as the sum of variables with increasing integer power and their coefficients in a certain field. For example the following might be still known from school:

```
P(x) = x^2 + 4x + 3
```

## Examples


Adding two polynomials
---
```javascript
const p = new Polynomial("3x^2").add("-x^2"); // 2x^2
```

Second derivative of polynomial
---
```javascript
const p = new Polynomial("5+3x^3+6x^5").derive(2); // 120x^3+18x
```

## Parser


Any function (see below) as well as the constructor of the *Polynomial* class parses its input like this.

You can pass either Objects, Doubles or Strings. Make sure strings don't contain any white-spaces or brackets. The parser doesn't analyse the string recursively.

Objects
---
```javascript
new Polynomial({'3': 4, '5': '9'}); // 9x^5+4x^3
new Polynomial([1,2,3]); //3x^2+2x+1
```

Doubles
---
```javascript
new Polynomial(55); // 55x^0
```

Strings
---
```javascript
new Polynomial("98x^2+4+23x^4");
```
The string parser passes every coefficient directly to the field parser, which allows to pass complex and rational coefficients as well:

```javascript
// Example with complex numbers
Polynomial.setField("C");
new Polynomial("98x^2+i+23ix^4");

// Example with rational numbers
Polynomial.setField("Q");
new Polynomial("5/3x^3+4/3x");
```

## Fields


Polynomial.js is held general in order to operate on various fields. [Fraction.js](https://github.com/infusion/Fraction.js) and [Complex.js](https://github.com/infusion/Complex.js) build the perfect base to extend polynomials to rational and complex numbers.

* ℚ: Rational numbers supported by [Fraction.js](https://github.com/infusion/Fraction.js) 
* ℂ: Complex numbers supported by [Complex.js](https://github.com/infusion/Complex.js)
* ℍ: Quaternions supported by [Quaternion.js](https://github.com/infusion/Quaternion.js)
* ℤ<sub>p</sub>: Field of integers mod p, with p prime
* ℝ: Field of real numbers


## Accessing Coefficients

```javascript
const p = new Polynomial("98x^2+4+23x^4");

console.log(p.coeff);
```


### Examples

```javascript
Polynomial.setField("Q");
Polynomial("3/2x^2-4x").mod("5x"); // 0

Polynomial.setField("Z11");
Polynomial("3x^2+x+7").gcd("3x^2+x+7"); // x^2+4x+6

Polynomial.setField("Z7");
Polynomial("9x^2+4").pow(3); // x^6+6x^4+5x^2+1

Polynomial.setField("R");
Polynomial("3x^3-1").mul(4); // 12x^3-4

// Derivative of polynomial
Polynomial.setField("Q");
Polynomial("5+3x^3+6x^5").derive(); // 30x^4+9x^2

// Integrated polynomial
Polynomial.setField("Q");
Polynomial("3x^2").integrate(); // x^3
```

## Functions


Polynomial add(n)
---
Returns the sum of the actual polynomial and the parameter n

Polynomial sub(n)
---
Returns the difference of the actual polynomial and the parameter n

Polynomial mul(n)
---
Returns the product of the actual polynomial and the parameter n

Polynomial addmul(x, y)
---
Adds the product of x and y to the actual number

Polynomial div(n)
---
Returns the quotient of the actual polynomial and the parameter n

There is a global variable to enable division tracing like this, if you want to output details:

```javascript
Polynomial.trace = true;
new Polynomial("x^4+3x^3+2x^2+6x")
        .div("x+3");
console.log(Polynomial.trace.map(x => x.toString())); // ["x^4+3x^3", "2x^2+6x", "0"]
```

Polynomial neg(n)
---
Returns the negated polynomial

Polynomial reciprocal(n)
---
Returns the [reciprocal polynomial](https://en.wikipedia.org/wiki/Reciprocal_polynomial)

Polynomial lc()
---
Gets the leading coefficient

Polynomial lm()
---
Gets the leading monomial

Polynomial monic()
---
Divide all coefficients of the polynomial by lc()

Polynomial derive(n)
---
Returns the n-th derivative of the polynomial

Polynomial integrate(n)
---
Returns the n-th integration of the polynomial

mixed eval(x)
---
Evaluate the polynomial at point x, using [Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method). Type for x must be a valid value for the given field.

mixed result(x)
---
(Deprecated) Alias for `eval`.

Polynomial pow(exp)
---
Returns the power of the actual polynomial, raised to an integer exponent.

Polynomial mod(n)
---
Returns the modulus (rest of the division) of the actual polynomial and n (this % n).

Polynomial gcd(n)
---
Returns the greatest common divisor of two polynomials

Number degree()
---
Returns the degree of the polynomial

String toString()
---
Generates a string representation of the actual polynomial. This makes use of the `toString()` function of the field.

String toLatex()
---
Generates a LaTeX representation of the actual polynomial.

String toHorner()
---
Formats the actual polynomial to a [Horner Scheme](https://en.wikipedia.org/wiki/Horner's_method)

Polynomial clone()
---
Creates a copy of the actual Polynomial object

Polynomial Polynomial::fromRoots(roots)
---
Creates a new (monic) Polynomial whose roots lie at the values provided in the array `roots`

Polynomial::setField(x)
---
Sets the field globally. Choose one of the following strings for `x`:

- "R": real numbers
- "Q": rational numbers
- "H": quaternions
- "C": complex numbers
- "Zp": with p being a prime number, like "Z7"
- or an object with the field operators. See examples folders for BigInt

## Exceptions

If a really hard error occurs (parsing error, division by zero), *polynomial.js* throws exceptions! Please make sure you handle them correctly.


## Installation

Installing Polynomial.js is as easy as cloning this repo or use one of the following command:


```
npm install polynomial
```

## Using Polynomial.js with the browser

```html
<script src="fraction.min.js"></script> <!-- Needed for field/ring Q -->
<script src="complex.min.js"></script> <!-- Needed for field C -->
<script src="polynomial.min.js"></script>
<script>
Polynomial.setField("C")
console.log(Polynomial("4x+3i"));
</script>
```


## Coding Style

As every library I publish, Polynomial.js is also built to be as small as possible after compressing it with Google Closure Compiler in advanced mode. Thus the coding style orientates a little on maxing-out the compression rate. Please make sure you keep this style if you plan to extend the library.

## Building the library

After cloning the Git repository run:

```
npm install
npm run build
```

## Run a test

Testing the source against the shipped test suite is as easy as

```
npm run test
```

## Copyright and Licensing

Copyright (c) 2024, [Robert Eisele](https://raw.org/)
Licensed under the MIT license.
