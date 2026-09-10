'use strict';

const Polynomial = require('polynomial');

Polynomial.setField('Q');

// Newton iteration over rational coefficients. This example is experimental.
function toNumber(value) {
  const scale = 10n ** 15n;
  return Number(value.s * value.n * scale / value.d) / Number(scale);
}

function factor(input) {
  let polynomial = new Polynomial(input);
  let derivative = polynomial.derive();
  const roots = [];

  do {
    let previous;
    let root = new Polynomial(1).coeff[0];

    do {
      previous = root;
      root = root.sub(polynomial.eval(root).div(derivative.eval(root)));
    } while (Math.abs(toNumber(root) - toNumber(previous)) > 1e-10);

    const rootValue = toNumber(root);
    const nearestInteger = Math.round(rootValue);
    if (Math.abs(rootValue - nearestInteger) < 1e-7) {
      root = new Polynomial(nearestInteger).coeff[0];
    }

    polynomial = polynomial.div(new Polynomial('x').sub(root));
    derivative = polynomial.derive();
    roots.push(root);
  } while (derivative.degree() >= 0);

  return roots;
}

// Results are approximately 3, 4, and 9.
console.log(factor('x^3-16x^2+75x-108').map(root => root.toString()));
