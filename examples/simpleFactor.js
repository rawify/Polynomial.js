'use strict';

const Polynomial = require('polynomial');

Polynomial.setField('C');

const polynomial = new Polynomial('ix^2-13ix+36i');
const leadingCoefficient = polynomial.coeff[2];
const linearCoefficient = polynomial.coeff[1];
const constant = polynomial.coeff[0];
const discriminant = linearCoefficient
  .mul(linearCoefficient)
  .sub(leadingCoefficient.mul(constant).mul(4))
  .sqrt();
const denominator = leadingCoefficient.mul(2);
const firstRoot = linearCoefficient.neg().add(discriminant).div(denominator);
const secondRoot = linearCoefficient.neg().sub(discriminant).div(denominator);

// Outputs (X-9)(X-4).
console.log(`(X${firstRoot.neg()})(X${secondRoot.neg()})`);
