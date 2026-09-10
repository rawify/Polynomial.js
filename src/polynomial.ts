/**
 * @license Polynomial.js v2.0.0 9/10/2026
 * https://github.com/rawify/Polynomial.js
 *
 * Copyright (c) 2025-2026, Robert Eisele (https://raw.org/)
 * Licensed under the MIT license.
 */

import Complex from 'complex.js';
import Fraction from 'fraction.js';
import Quaternion from 'quaternion';

export type BuiltInField = 'R' | 'Q' | 'C' | 'H' | `Z${number}`;

export interface PolynomialField<C> {
  add(left: C, right: unknown): C;
  sub(left: C, right: unknown): C;
  mul(left: C, right: unknown): C;
  div(left: C, right: unknown): C;
  parse(value: unknown): C;
  empty(value: C | undefined): boolean;
  pow(base: C, exponent: number): C;
  equals(left: C, right: unknown): boolean;
  readonly characteristic?: number;
  magnitude?(value: C): number;
  sqrt?(value: C): C | null;
}

export type CoefficientMap<C> = Readonly<Record<number, C>>;

export type PolynomialInput<C> =
  | Polynomial<C>
  | CoefficientMap<C>
  | readonly unknown[]
  | C
  | string
  | number
  | null
  | undefined;

export type PolynomialResidueInput<C> = PolynomialResidue<C> | PolynomialInput<C>;

export interface PolynomialQuotientRing<C> {
  readonly modulus: Polynomial<C>;
  create(input?: PolynomialResidueInput<C>): PolynomialResidue<C>;
  zero(): PolynomialResidue<C>;
  one(): PolynomialResidue<C>;
  evaluate(
    polynomial: PolynomialInput<C>,
    value: PolynomialResidueInput<C>,
  ): PolynomialResidue<C>;
  isRoot(
    polynomial: PolynomialInput<C>,
    value: PolynomialResidueInput<C>,
  ): boolean;
  order(): number;
  isField(): boolean;
  elements(maxElements?: number): PolynomialResidue<C>[];
  primitiveElements(maxElements?: number): PolynomialResidue<C>[];
  primitivePolynomials(maxElements?: number): Polynomial<C>[];
}

export interface PolynomialResidue<C> {
  readonly ring: PolynomialQuotientRing<C>;
  readonly value: Polynomial<C>;
  add(input: PolynomialResidueInput<C>): PolynomialResidue<C>;
  sub(input: PolynomialResidueInput<C>): PolynomialResidue<C>;
  mul(input: PolynomialResidueInput<C>): PolynomialResidue<C>;
  div(input: PolynomialResidueInput<C>): PolynomialResidue<C>;
  neg(): PolynomialResidue<C>;
  pow(exponent: number): PolynomialResidue<C>;
  inverse(): PolynomialResidue<C>;
  additiveOrder(): number;
  additiveSubgroup(maxElements?: number): PolynomialResidue<C>[];
  multiplicativeOrder(): number;
  multiplicativeSubgroup(maxElements?: number): PolynomialResidue<C>[];
  isPrimitive(): boolean;
  minimalPolynomial(): Polynomial<C>;
  equals(input: PolynomialResidueInput<C>): boolean;
  isZero(): boolean;
  toString(): string;
}

type MutableCoefficientMap<C> = Record<number, C>;
type AlgebraicValue = {
  add(value: unknown): AlgebraicValue;
  sub(value: unknown): AlgebraicValue;
  mul(value: unknown): AlgebraicValue;
  div(value: unknown): AlgebraicValue;
  pow(value: unknown): AlgebraicValue;
  equals(value: unknown): boolean;
  sqrt?: () => AlgebraicValue;
  valueOf(): unknown;
  toString(): string;
  toLatex?: () => string;
};
type AlgebraicConstructor = new (value?: unknown) => AlgebraicValue;

const TERM_PATTERN = /([+-]?)(?:([^+x-]+)?(?:x(?:\^(\d+))?)|([^+x-]+))/g;

const realField: PolynomialField<number> = {
  add: (left, right) => left + Number(right),
  sub: (left, right) => left - Number(right),
  mul: (left, right) => left * Number(right),
  div: (left, right) => {
    const divisor = Number(right);
    if (divisor === 0) {
      throw new RangeError('Division by zero');
    }
    return left / divisor;
  },
  parse: value => Number.parseFloat(String(value)),
  empty: value => value === undefined || value === 0,
  pow: (base, exponent) => base ** exponent,
  equals: (left, right) => left === Number(right),
  characteristic: 0,
  magnitude: value => Math.abs(value),
  sqrt: value => {
    const root = Math.sqrt(value);
    return Number.isNaN(root) ? null : root;
  },
};

let activeField: PolynomialField<unknown> = realField as PolynomialField<unknown>;

function createCoefficientMap<C>(): MutableCoefficientMap<C> {
  return Object.create(null) as MutableCoefficientMap<C>;
}

function cloneCoefficientMap<C>(source: CoefficientMap<C>): MutableCoefficientMap<C> {
  return Object.assign(createCoefficientMap<C>(), source);
}

function normalize<C>(coefficients: MutableCoefficientMap<C>, field: PolynomialField<C>): MutableCoefficientMap<C> {
  for (const exponent of Object.keys(coefficients)) {
    const numericExponent = Number(exponent);
    if (field.empty(coefficients[numericExponent])) {
      delete coefficients[numericExponent];
    }
  }
  return coefficients;
}

function degreeOf<C>(coefficients: CoefficientMap<C>, field: PolynomialField<C>): number {
  let degree = -Infinity;
  for (const exponent of Object.keys(coefficients)) {
    const numericExponent = Number(exponent);
    if (!field.empty(coefficients[numericExponent])) {
      degree = Math.max(degree, numericExponent);
    }
  }
  return degree;
}

function leadingExponent<C>(coefficients: CoefficientMap<C>, field: PolynomialField<C>): number | null {
  const degree = degreeOf(coefficients, field);
  return Number.isFinite(degree) ? degree : null;
}

function makeMonic<C>(coefficients: MutableCoefficientMap<C>, field: PolynomialField<C>): MutableCoefficientMap<C> {
  const exponent = leadingExponent(coefficients, field);
  if (exponent === null) {
    return coefficients;
  }

  const leadingCoefficient = coefficients[exponent];
  for (const key of Object.keys(coefficients)) {
    const numericExponent = Number(key);
    coefficients[numericExponent] = field.div(
      coefficients[numericExponent] as C,
      leadingCoefficient,
    );
  }
  return normalize(coefficients, field);
}

function modulo(value: number, modulus: number): number {
  return ((value % modulus) + modulus) % modulus;
}

function isPrime(value: number): boolean {
  if (!Number.isSafeInteger(value) || value < 2) {
    return false;
  }
  const candidate = BigInt(value);
  for (const prime of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n]) {
    if (candidate === prime) {
      return true;
    }
    if (candidate % prime === 0n) {
      return false;
    }
  }

  let oddPart = candidate - 1n;
  let powerOfTwo = 0;
  while (oddPart % 2n === 0n) {
    oddPart /= 2n;
    powerOfTwo += 1;
  }

  const power = (base: bigint, exponent: bigint): bigint => {
    let result = 1n;
    let factor = base % candidate;
    let remaining = exponent;
    while (remaining > 0n) {
      if (remaining % 2n === 1n) {
        result = result * factor % candidate;
      }
      factor = factor * factor % candidate;
      remaining /= 2n;
    }
    return result;
  };

  return [2n, 325n, 9375n, 28178n, 450775n, 9780504n, 1795265022n]
    .every(base => {
      if (base % candidate === 0n) {
        return true;
      }
      let witness = power(base, oddPart);
      if (witness === 1n || witness === candidate - 1n) {
        return true;
      }
      for (let index = 1; index < powerOfTwo; index += 1) {
        witness = witness * witness % candidate;
        if (witness === candidate - 1n) {
          return true;
        }
      }
      return false;
    });
}

function primeDivisors(value: number): number[] {
  const divisors: number[] = [];
  let remainder = value;
  for (let divisor = 2; divisor * divisor <= remainder; divisor += divisor === 2 ? 1 : 2) {
    if (remainder % divisor !== 0) {
      continue;
    }
    divisors.push(divisor);
    while (remainder % divisor === 0) {
      remainder /= divisor;
    }
  }
  if (remainder > 1) {
    divisors.push(remainder);
  }
  return divisors;
}

function extendedGcd(left: number, right: number): [number, number, number] {
  let oldRemainder = left;
  let remainder = right;
  let oldCoefficient = 1;
  let coefficient = 0;
  let oldOtherCoefficient = 0;
  let otherCoefficient = 1;

  while (remainder !== 0) {
    const quotient = Math.trunc(oldRemainder / remainder);
    [oldRemainder, remainder] = [remainder, oldRemainder - quotient * remainder];
    [oldCoefficient, coefficient] = [coefficient, oldCoefficient - quotient * coefficient];
    [oldOtherCoefficient, otherCoefficient] = [
      otherCoefficient,
      oldOtherCoefficient - quotient * otherCoefficient,
    ];
  }

  return [Math.abs(oldRemainder), oldCoefficient, oldOtherCoefficient];
}

function modularInverse(value: number, modulus: number): number {
  const normalized = modulo(value, modulus);
  const [divisor, coefficient] = extendedGcd(normalized, modulus);
  if (divisor !== 1) {
    throw new RangeError(`${value} has no multiplicative inverse modulo ${modulus}`);
  }
  return modulo(coefficient, modulus);
}

function modularPower(base: number, exponent: number, modulus: number): number {
  let result = 1;
  let factor = modulo(base, modulus);
  let power = exponent;

  while (power > 0) {
    if (power % 2 === 1) {
      result = modulo(result * factor, modulus);
    }
    factor = modulo(factor * factor, modulus);
    power = Math.floor(power / 2);
  }
  return result;
}

function modularSquareRoot(value: number, modulus: number): number | null {
  const normalized = modulo(value, modulus);
  if (normalized === 0 || modulus === 2) {
    return normalized;
  }
  if (modulus % 2 === 0) {
    return null;
  }
  if (modularPower(normalized, (modulus - 1) / 2, modulus) !== 1) {
    return null;
  }
  if (modulus % 4 === 3) {
    const root = modularPower(normalized, (modulus + 1) / 4, modulus);
    return Math.min(root, modulus - root);
  }

  let oddPart = modulus - 1;
  let powerOfTwo = 0;
  while (oddPart % 2 === 0) {
    oddPart /= 2;
    powerOfTwo += 1;
  }

  let nonResidue = 2;
  while (
    nonResidue < modulus
    && modularPower(nonResidue, (modulus - 1) / 2, modulus) !== modulus - 1
  ) {
    nonResidue += 1;
  }
  if (nonResidue === modulus) {
    return null;
  }

  let order = powerOfTwo;
  let factor = modularPower(nonResidue, oddPart, modulus);
  let remainder = modularPower(normalized, oddPart, modulus);
  let root = modularPower(normalized, (oddPart + 1) / 2, modulus);

  while (remainder !== 1) {
    let index = 1;
    let square = modulo(remainder * remainder, modulus);
    while (index < order && square !== 1) {
      square = modulo(square * square, modulus);
      index += 1;
    }
    if (index === order) {
      return null;
    }

    const correction = modularPower(factor, 2 ** (order - index - 1), modulus);
    root = modulo(root * correction, modulus);
    factor = modulo(correction * correction, modulus);
    remainder = modulo(remainder * factor, modulus);
    order = index;
  }

  return Math.min(root, modulus - root);
}

function createModularField(modulus: number): PolynomialField<number> {
  if (!Number.isSafeInteger(modulus) || modulus <= 1) {
    throw new RangeError('A modular field requires an integer modulus greater than one');
  }

  return {
    add: (left, right) => modulo(left + Number(right), modulus),
    sub: (left, right) => modulo(left - Number(right), modulus),
    mul: (left, right) => modulo(left * Number(right), modulus),
    div: (left, right) => modulo(left * modularInverse(Number(right), modulus), modulus),
    parse: value => modulo(Number.parseInt(String(value), 10), modulus),
    empty: value => value === undefined || value === 0,
    pow: (base, exponent) => modularPower(base, exponent, modulus),
    equals: (left, right) => left === modulo(Number(right), modulus),
    characteristic: modulus,
    sqrt: value => modularSquareRoot(value, modulus),
  };
}

function createAlgebraicField(constructor: AlgebraicConstructor): PolynomialField<AlgebraicValue> {
  const field: PolynomialField<AlgebraicValue> = {
    add: (left, right) => new constructor(left).add(right),
    sub: (left, right) => new constructor(left).sub(right),
    mul: (left, right) => new constructor(left).mul(right),
    div: (left, right) => new constructor(left).div(right),
    parse: value => new constructor(value),
    empty: value => value === undefined || new constructor(value).equals(0),
    pow: (base, exponent) => new constructor(base).pow(exponent),
    equals: (left, right) => new constructor(left).equals(right),
    characteristic: 0,
  };

  const FractionClass = Fraction as unknown as AlgebraicConstructor;
  if (constructor === FractionClass) {
    field.magnitude = value => {
      const fraction = new constructor(value) as AlgebraicValue & { n: bigint; d: bigint };
      return bigintRatioMagnitude(fraction.n, fraction.d);
    };
    field.sqrt = value => {
      const fraction = new constructor(value) as AlgebraicValue & {
        s: bigint;
        n: bigint;
        d: bigint;
      };
      if (fraction.s < 0n) {
        return null;
      }

      const numerator = bigintSquareRoot(fraction.n);
      const denominator = bigintSquareRoot(fraction.d);
      if (numerator * numerator !== fraction.n || denominator * denominator !== fraction.d) {
        return null;
      }
      return new constructor({ s: 1n, n: numerator, d: denominator });
    };
  } else if (constructor === (Complex as unknown as AlgebraicConstructor)) {
    field.magnitude = value => {
      const complex = new constructor(value) as AlgebraicValue & { abs(): number };
      return complex.abs();
    };
    field.sqrt = value => new constructor(value).sqrt?.() ?? null;
  }

  return field;
}

function bigintSquareRoot(value: bigint): bigint {
  if (value < 2n) {
    return value;
  }

  let estimate = 1n << (BigInt(value.toString(2).length) + 1n) / 2n;
  let next = (estimate + value / estimate) / 2n;
  while (next < estimate) {
    estimate = next;
    next = (estimate + value / estimate) / 2n;
  }
  return estimate;
}

function bigintRatioMagnitude(numerator: bigint, denominator: bigint): number {
  if (numerator === 0n) {
    return 0;
  }

  const numeratorDigits = (numerator < 0n ? -numerator : numerator).toString();
  const denominatorDigits = (denominator < 0n ? -denominator : denominator).toString();
  const precision = 16;
  const numeratorPrefix = numeratorDigits.slice(0, precision);
  const denominatorPrefix = denominatorDigits.slice(0, precision);
  const numeratorMantissa = Number(numeratorPrefix) / 10 ** (numeratorPrefix.length - 1);
  const denominatorMantissa = Number(denominatorPrefix) / 10 ** (denominatorPrefix.length - 1);

  return numeratorMantissa
    / denominatorMantissa
    * 10 ** (numeratorDigits.length - denominatorDigits.length);
}

function fractionToNumber(fraction: { s: bigint; n: bigint; d: bigint }): number {
  const numerator = Number(fraction.n);
  const denominator = Number(fraction.d);
  if (Number.isFinite(numerator) && Number.isFinite(denominator)) {
    return Number(fraction.s) * numerator / denominator;
  }
  return Number(fraction.s) * bigintRatioMagnitude(fraction.n, fraction.d);
}

function binomialCoefficient(total: number, selected: number): number {
  const count = Math.min(selected, total - selected);
  let result = 1;
  for (let index = 1; index <= count; index += 1) {
    result = result * (total - count + index) / index;
  }
  return result;
}

function assertField<C>(field: PolynomialField<C>): void {
  const methods: Array<keyof PolynomialField<C>> = [
    'add',
    'sub',
    'mul',
    'div',
    'parse',
    'empty',
    'pow',
    'equals',
  ];

  for (const method of methods) {
    if (typeof field[method] !== 'function') {
      throw new TypeError(`Polynomial field is missing ${method}()`);
    }
  }
}

function parseExpression<C>(sign: string, expression: string, field: PolynomialField<C>): C {
  const tokens = expression.match(/[^*/]+|[*/]/g);
  if (tokens === null || tokens.length === 0) {
    throw new SyntaxError('Invalid coefficient expression');
  }

  let result = field.parse(`${sign}${tokens[0]}`);
  for (let index = 1; index < tokens.length; index += 2) {
    const operator = tokens[index];
    const operand = field.parse(tokens[index + 1] ?? 1);
    if (operator === '*') {
      result = field.mul(result, operand);
    } else if (operator === '/') {
      result = field.div(result, operand);
    } else {
      throw new SyntaxError(`Invalid coefficient operator: ${operator}`);
    }
  }
  return result;
}

function isAlgebraicConstant(value: unknown): boolean {
  const FractionClass = Fraction as unknown as AlgebraicConstructor;
  const QuaternionClass = Quaternion as unknown as AlgebraicConstructor;
  return value instanceof FractionClass
    || value instanceof Complex
    || value instanceof QuaternionClass;
}

function parsePolynomial<C>(input: PolynomialInput<C>, field: PolynomialField<C>): MutableCoefficientMap<C> {
  const coefficients = createCoefficientMap<C>();
  if (input === null || input === undefined) {
    return coefficients;
  }

  if (typeof input === 'number') {
    const coefficient = field.parse(input);
    if (!field.empty(coefficient)) {
      coefficients[0] = coefficient;
    }
    return coefficients;
  }

  if (typeof input === 'string') {
    const source = input.replace(/\s+/g, '');
    if (source.length === 0) {
      return coefficients;
    }

    TERM_PATTERN.lastIndex = 0;
    let cursor = 0;
    let match: RegExpExecArray | null;
    while ((match = TERM_PATTERN.exec(source)) !== null) {
      if (match.index !== cursor || match[0].length === 0) {
        throw new SyntaxError(`Invalid polynomial near "${source.slice(cursor)}"`);
      }
      if (cursor > 0 && match[1] === '') {
        throw new SyntaxError(`Missing operator near "${source.slice(cursor)}"`);
      }

      const constant = match[4];
      const coefficientExpression = constant ?? match[2] ?? '1';
      const exponent = constant === undefined ? Number.parseInt(match[3] ?? '1', 10) : 0;
      const coefficient = parseExpression(match[1] ?? '', coefficientExpression, field);
      const previous = coefficients[exponent];
      coefficients[exponent] = previous === undefined
        ? coefficient
        : field.add(previous, coefficient);
      cursor = TERM_PATTERN.lastIndex;
    }

    if (cursor !== source.length) {
      throw new SyntaxError(`Invalid polynomial near "${source.slice(cursor)}"`);
    }
    return normalize(coefficients, field);
  }

  if (typeof input === 'object') {
    if (isAlgebraicConstant(input)) {
      const coefficient = field.parse(input);
      if (!field.empty(coefficient)) {
        coefficients[0] = coefficient;
      }
      return coefficients;
    }

    const source = input instanceof Polynomial
      ? input.coeff
      : input;
    for (const [key, value] of Object.entries(source)) {
      const exponent = Number(key);
      if (!Number.isSafeInteger(exponent) || exponent < 0) {
        throw new RangeError(`Invalid polynomial exponent: ${key}`);
      }

      const coefficient = field.parse(value);
      if (!field.empty(coefficient)) {
        coefficients[exponent] = coefficient;
      }
    }
    return coefficients;
  }

  throw new TypeError('Invalid polynomial input');
}

function fallingFactorial(value: number, count: number): number {
  let result = 1;
  for (let index = 0; index < count; index += 1) {
    result *= value - index;
  }
  return result;
}

function risingFactorial(value: number, count: number): number {
  let result = 1;
  for (let index = 1; index <= count; index += 1) {
    result *= value + index;
  }
  return result;
}

function isComplex(value: unknown): value is Complex {
  return value instanceof Complex;
}

export class Polynomial<C = number> {
  static trace: boolean | Polynomial<unknown>[] | null = null;

  readonly coeff: MutableCoefficientMap<C>;
  private readonly field: PolynomialField<C>;

  constructor(
    input?: PolynomialInput<C>,
    field: PolynomialField<C> = activeField as PolynomialField<C>,
  ) {
    assertField(field);
    this.field = field;
    this.coeff = parsePolynomial(input, field);
  }

  gcd(input: PolynomialInput<C>): Polynomial<C> {
    let left = normalize(cloneCoefficientMap(this.coeff), this.field);
    let right = parsePolynomial(input, this.field);

    if (degreeOf(left, this.field) < 0) {
      return this.create(makeMonic(right, this.field));
    }
    if (degreeOf(right, this.field) < 0) {
      return this.create(makeMonic(left, this.field));
    }

    while (degreeOf(right, this.field) >= 0) {
      const remainder = cloneCoefficientMap(left);
      this.divideCoefficients(remainder, right);
      left = right;
      right = normalize(remainder, this.field);
    }

    return this.create(makeMonic(left, this.field));
  }

  lcm(input: PolynomialInput<C>): Polynomial<C> {
    const right = this.create(parsePolynomial(input, this.field));
    if (this.isZero() || right.isZero()) {
      return this.create({});
    }
    return this.mul(right).div(this.gcd(right)).monic();
  }

  neg(): Polynomial<C> {
    const result = createCoefficientMap<C>();
    for (const [key, coefficient] of Object.entries(this.coeff)) {
      result[Number(key)] = this.field.mul(coefficient, -1);
    }
    return this.create(result);
  }

  reciprocal(): Polynomial<C> {
    const result = createCoefficientMap<C>();
    const degree = this.degree();
    if (degree < 0) {
      return this.create(result);
    }

    for (const [key, coefficient] of Object.entries(this.coeff)) {
      result[degree - Number(key)] = coefficient;
    }
    return this.create(result);
  }

  sqrt(): Polynomial<C> | null {
    const degree = this.degree();
    if (degree < 0) {
      return this.create({});
    }
    if (degree % 2 !== 0) {
      return null;
    }
    if (this.field === (realField as unknown as PolynomialField<C>)) {
      const exactCoefficients = createCoefficientMap<AlgebraicValue>();
      for (const [key, coefficient] of Object.entries(this.coeff)) {
        exactCoefficients[Number(key)] = new (Fraction as unknown as AlgebraicConstructor)(
          String(coefficient),
        );
      }
      const rationalField = createAlgebraicField(
        Fraction as unknown as AlgebraicConstructor,
      );
      const exactRoot = new Polynomial<AlgebraicValue>(exactCoefficients, rationalField).sqrt();
      if (exactRoot !== null) {
        const coefficients = createCoefficientMap<C>();
        for (const [key, coefficient] of Object.entries(exactRoot.coeff)) {
          const fraction = coefficient as AlgebraicValue & { s: bigint; n: bigint; d: bigint };
          coefficients[Number(key)] = this.field.parse(fractionToNumber(fraction));
        }
        return this.create(coefficients);
      }
    }
    if (this.field.sqrt === undefined) {
      throw new TypeError('Polynomial square root requires a coefficient square root function');
    }

    if (this.field.characteristic === 2) {
      const coefficients = createCoefficientMap<C>();
      for (const [key, coefficient] of Object.entries(this.coeff)) {
        const exponent = Number(key);
        if (exponent % 2 !== 0) {
          return null;
        }
        const root = this.field.sqrt(coefficient);
        if (root === null) {
          return null;
        }
        if (!this.field.empty(root)) {
          coefficients[exponent / 2] = root;
        }
      }
      const candidate = this.create(coefficients);
      return this.isSquareRoot(candidate) ? candidate : null;
    }

    const rootDegree = degree / 2;
    const leadingRoot = this.field.sqrt(this.lc() as C);
    if (leadingRoot === null || this.field.empty(leadingRoot)) {
      return null;
    }

    const coefficients = createCoefficientMap<C>();
    coefficients[rootDegree] = leadingRoot;
    for (let exponent = rootDegree - 1; exponent >= 0; exponent -= 1) {
      const sourceExponent = rootDegree + exponent;
      let remainder = this.coeff[sourceExponent] ?? this.field.parse(0);

      for (let leftExponent = exponent + 1; leftExponent < rootDegree; leftExponent += 1) {
        const rightExponent = sourceExponent - leftExponent;
        const left = coefficients[leftExponent];
        const right = coefficients[rightExponent];
        if (left !== undefined && right !== undefined) {
          remainder = this.field.sub(remainder, this.field.mul(left, right));
        }
      }

      const coefficient = this.field.div(
        remainder,
        this.field.mul(leadingRoot, 2),
      );
      if (!this.field.empty(coefficient)) {
        coefficients[exponent] = coefficient;
      }
    }

    const candidate = this.create(coefficients);
    return this.isSquareRoot(candidate) ? candidate : null;
  }

  eval(value: unknown): C | number {
    const degree = this.degree();
    if (degree < 0) {
      return 0;
    }

    let result = this.coeff[degree] as C;
    for (let exponent = degree - 1; exponent >= 0; exponent -= 1) {
      result = this.field.mul(result, value);
      const coefficient = this.coeff[exponent];
      if (!this.field.empty(coefficient)) {
        result = this.field.add(result, coefficient);
      }
    }
    return result;
  }

  compose(input: PolynomialInput<C>): Polynomial<C> {
    const inner = this.create(parsePolynomial(input, this.field));
    let result = this.create({});

    for (let exponent = this.degree(); exponent >= 0; exponent -= 1) {
      result = result.mul(inner);
      const coefficient = this.coeff[exponent];
      if (!this.field.empty(coefficient)) {
        result = result.add({ 0: coefficient as C });
      }
    }
    return result;
  }

  /** @deprecated Use eval(). */
  result(value: unknown): C | number {
    return this.eval(value);
  }

  lc(): C | undefined {
    const exponent = leadingExponent(this.coeff, this.field);
    return exponent === null ? undefined : this.coeff[exponent];
  }

  lm(): Polynomial<C> {
    const result = createCoefficientMap<C>();
    const exponent = leadingExponent(this.coeff, this.field);
    if (exponent !== null) {
      result[exponent] = this.coeff[exponent] as C;
    }
    return this.create(result);
  }

  monic(): Polynomial<C> {
    return this.create(makeMonic(cloneCoefficientMap(this.coeff), this.field));
  }

  equals(input: PolynomialInput<C>): boolean {
    const right = parsePolynomial(input, this.field);
    const degree = Math.max(this.degree(), degreeOf(right, this.field));
    const zero = this.field.parse(0);

    for (let exponent = 0; exponent <= degree; exponent += 1) {
      if (!this.field.equals(this.coeff[exponent] ?? zero, right[exponent] ?? zero)) {
        return false;
      }
    }
    return true;
  }

  isZero(): boolean {
    return this.degree() < 0;
  }

  isMonic(): boolean {
    const leadingCoefficient = this.lc();
    return leadingCoefficient !== undefined && this.field.equals(leadingCoefficient, 1);
  }

  compare(input: PolynomialInput<C>): -1 | 0 | 1 {
    this.assertPrimeField();
    const right = this.create(parsePolynomial(input, this.field));
    const degree = Math.max(this.degree(), right.degree());

    for (let exponent = degree; exponent >= 0; exponent -= 1) {
      const difference = this.coefficientDigit(exponent) - right.coefficientDigit(exponent);
      if (difference !== 0) {
        return difference < 0 ? -1 : 1;
      }
    }
    return 0;
  }

  successor(): Polynomial<C> {
    this.assertPrimeField();
    const result = cloneCoefficientMap(this.coeff);
    let exponent = 0;

    while (true) {
      const coefficient = result[exponent] ?? this.field.parse(0);
      const next = this.field.add(coefficient, 1);
      if (!this.field.empty(next)) {
        result[exponent] = next;
        break;
      }
      delete result[exponent];
      exponent += 1;
    }
    return this.create(result);
  }

  isDivisibleBy(input: PolynomialInput<C>): boolean {
    return this.mod(input).isZero();
  }

  factor(): Polynomial<C>[] {
    this.assertPrimeField();
    if (this.isZero()) {
      throw new RangeError('The zero polynomial does not have a finite factorization');
    }

    const leadingCoefficient = this.lc() as C;
    const factors = this.degree() > 0
      ? this.factorMonicFinite(this.monic())
      : [];
    factors.sort((left, right) => (
      left.degree() - right.degree()
      || left.toString().localeCompare(right.toString())
    ));

    if (!this.field.equals(leadingCoefficient, 1)) {
      factors.unshift(this.create({ 0: leadingCoefficient }));
    }
    return factors;
  }

  isIrreducible(): boolean {
    this.assertPrimeField();
    if (this.degree() < 1) {
      return false;
    }
    const factors = this.factorMonicFinite(this.monic());
    return factors.length === 1 && factors[0]?.degree() === this.degree();
  }

  cyclicCodeBasis(length: number): Polynomial<C>[] {
    this.assertPrimeField();
    if (!Number.isSafeInteger(length) || length <= 0) {
      throw new RangeError('Cyclic code length must be a positive safe integer');
    }
    if (this.isZero()) {
      throw new RangeError('Cyclic code generator must be non-zero');
    }

    const generator = this.monic();
    if (generator.degree() > length) {
      throw new RangeError('Cyclic code generator degree exceeds the code length');
    }
    const ambient = this.create({
      0: this.field.mul(this.field.parse(1), -1),
      [length]: this.field.parse(1),
    });
    if (!ambient.isDivisibleBy(generator)) {
      throw new RangeError('Cyclic code generator must divide x^length-1');
    }

    return Array.from(
      { length: length - generator.degree() },
      (_, shift) => generator.mul({ [shift]: this.field.parse(1) }),
    );
  }

  dispersion(maxShift?: number): number {
    if (maxShift !== undefined && (!Number.isSafeInteger(maxShift) || maxShift < 0)) {
      throw new RangeError('Maximum dispersion shift must be a non-negative safe integer');
    }
    if (this.degree() <= 0) {
      return 0;
    }
    if (this.field === (realField as unknown as PolynomialField<C>)) {
      const exactCoefficients = createCoefficientMap<AlgebraicValue>();
      for (const [key, coefficient] of Object.entries(this.coeff)) {
        exactCoefficients[Number(key)] = new (Fraction as unknown as AlgebraicConstructor)(
          String(coefficient),
        );
      }
      const rationalField = createAlgebraicField(
        Fraction as unknown as AlgebraicConstructor,
      );
      return new Polynomial<AlgebraicValue>(exactCoefficients, rationalField)
        .dispersion(maxShift);
    }

    let upperBound = maxShift;
    if (upperBound === undefined) {
      if ((this.field.characteristic ?? 0) > 0) {
        return Infinity;
      }
      if (this.field.magnitude === undefined) {
        throw new TypeError('Automatic dispersion requires a coefficient magnitude function');
      }

      const leadingCoefficient = this.lc() as C;
      let largestRatio = 0;
      for (const [key, coefficient] of Object.entries(this.coeff)) {
        if (Number(key) < this.degree()) {
          const ratio = this.field.div(coefficient, leadingCoefficient);
          largestRatio = Math.max(largestRatio, this.field.magnitude(ratio));
        }
      }
      if (!Number.isFinite(largestRatio)) {
        throw new RangeError('Unable to determine a finite dispersion bound');
      }

      upperBound = Math.ceil(2 * (1 + largestRatio) * (1 + 1e-12)) + 1;
      if (!Number.isSafeInteger(upperBound)) {
        throw new RangeError('Dispersion bound exceeds the safe integer range');
      }
    }

    for (let shift = upperBound; shift > 0; shift -= 1) {
      if (this.gcd(this.shift(shift)).degree() > 0) {
        return shift;
      }
    }
    return 0;
  }

  add(input: PolynomialInput<C>): Polynomial<C> {
    const right = parsePolynomial(input, this.field);
    const result = cloneCoefficientMap(this.coeff);

    for (const [key, coefficient] of Object.entries(right)) {
      const exponent = Number(key);
      const current = result[exponent];
      result[exponent] = current === undefined
        ? coefficient
        : this.field.add(current, coefficient);
    }
    return this.create(normalize(result, this.field));
  }

  sub(input: PolynomialInput<C>): Polynomial<C> {
    const right = parsePolynomial(input, this.field);
    const result = cloneCoefficientMap(this.coeff);

    for (const [key, coefficient] of Object.entries(right)) {
      const exponent = Number(key);
      const current = result[exponent];
      result[exponent] = current === undefined
        ? this.field.mul(coefficient, -1)
        : this.field.sub(current, coefficient);
    }
    return this.create(normalize(result, this.field));
  }

  mul(input: PolynomialInput<C>): Polynomial<C> {
    const right = parsePolynomial(input, this.field);
    const result = createCoefficientMap<C>();

    for (const [leftKey, leftCoefficient] of Object.entries(this.coeff)) {
      for (const [rightKey, rightCoefficient] of Object.entries(right)) {
        const exponent = Number(leftKey) + Number(rightKey);
        const product = this.field.mul(leftCoefficient, rightCoefficient);
        const current = result[exponent];
        result[exponent] = current === undefined
          ? product
          : this.field.add(current, product);
      }
    }
    return this.create(normalize(result, this.field));
  }

  addmul(left: PolynomialInput<C>, right: PolynomialInput<C>): Polynomial<C> {
    return this.add(this.create(parsePolynomial(left, this.field)).mul(right));
  }

  div(input: PolynomialInput<C>): Polynomial<C> {
    const numerator = cloneCoefficientMap(this.coeff);
    const denominator = parsePolynomial(input, this.field);
    return this.create(this.divideCoefficients(numerator, denominator));
  }

  mod(input: PolynomialInput<C>): Polynomial<C> {
    const remainder = cloneCoefficientMap(this.coeff);
    this.divideCoefficients(remainder, parsePolynomial(input, this.field));
    return this.create(normalize(remainder, this.field));
  }

  pow(exponent: number): Polynomial<C> {
    if (!Number.isSafeInteger(exponent) || exponent < 0) {
      throw new RangeError('Polynomial exponent must be a non-negative safe integer');
    }

    let result = this.create({ 0: this.field.parse(1) });
    let factor = this.clone();
    let power = exponent;
    while (power > 0) {
      if (power % 2 === 1) {
        result = result.mul(factor);
      }
      power = Math.floor(power / 2);
      if (power > 0) {
        factor = factor.mul(factor);
      }
    }
    return result;
  }

  derive(): Polynomial<C>;
  derive(order: number): Polynomial<C> | null;
  derive(order = 1): Polynomial<C> | null {
    if (!Number.isSafeInteger(order) || order < 0) {
      return null;
    }

    const result = createCoefficientMap<C>();
    for (const [key, coefficient] of Object.entries(this.coeff)) {
      const exponent = Number(key);
      if (exponent >= order) {
        result[exponent - order] = this.field.mul(
          coefficient,
          fallingFactorial(exponent, order),
        );
      }
    }
    return this.create(normalize(result, this.field));
  }

  integrate(): Polynomial<C>;
  integrate(order: number): Polynomial<C> | null;
  integrate(order = 1): Polynomial<C> | null {
    if (!Number.isSafeInteger(order) || order < 0) {
      return null;
    }

    const result = createCoefficientMap<C>();
    for (const [key, coefficient] of Object.entries(this.coeff)) {
      const exponent = Number(key);
      result[exponent + order] = this.field.div(
        coefficient,
        risingFactorial(exponent, order),
      );
    }
    return this.create(normalize(result, this.field));
  }

  degree(): number {
    return degreeOf(this.coeff, this.field);
  }

  toArray(width?: number): C[] {
    if (width !== undefined && (!Number.isSafeInteger(width) || width < 0)) {
      throw new RangeError('Coefficient array width must be a non-negative safe integer');
    }

    const degree = this.degree();
    const length = Math.max(width ?? 0, Number.isFinite(degree) ? degree + 1 : 0);
    return Array.from(
      { length },
      (_, exponent) => this.coeff[exponent] ?? this.field.parse(0),
    );
  }

  clone(): Polynomial<C> {
    return this.create(this.coeff);
  }

  toString(): string {
    return this.format('toString');
  }

  toLatex(): string {
    return this.format('toLatex');
  }

  toHorner(): string {
    const exponents = Object.keys(this.coeff)
      .map(Number)
      .filter(exponent => !this.field.empty(this.coeff[exponent]))
      .sort((left, right) => left - right);

    if (exponents.length === 0) {
      return '0';
    }

    const formatCoefficient = (coefficient: C, hasSign: boolean): string => {
      if (isComplex(coefficient)) {
        if (coefficient.im === 0) {
          return `${coefficient.re > 0 && hasSign ? '+' : ''}${coefficient.re}`;
        }
        if (coefficient.re === 0) {
          let result = '';
          if (coefficient.im === -1) {
            result = '-';
          } else if (coefficient.im !== 1) {
            result = `${coefficient.im > 0 && hasSign ? '+' : ''}${coefficient.im}`;
          } else if (hasSign) {
            result = '+';
          }
          return `${result}i`;
        }
        return `${hasSign ? '+' : ''}(${coefficient.toString()})`;
      }

      const text = String(coefficient);
      const numeric = Number((coefficient as { valueOf?: () => unknown })?.valueOf?.() ?? coefficient);
      return `${numeric > 0 && hasSign ? '+' : ''}${text}`;
    };

    const build = (position: number): string => {
      const index = exponents.length - position - 1;
      const exponent = exponents[index] as number;
      const previousExponent = exponents[index - 1] ?? 0;
      const difference = exponent - previousExponent;
      let factor = difference > 0 ? `x${difference > 1 ? `^${difference}` : ''}` : '';

      if (index > 0) {
        factor += formatCoefficient(this.coeff[previousExponent] as C, true);
      }
      if (position === 0) {
        return `${formatCoefficient(this.coeff[exponent] as C, false)}${factor}`;
      }

      const nested = build(position - 1);
      return `${exponent !== 0 ? `(${nested})` : nested}${factor}`;
    };

    return build(exponents.length - 1);
  }

  static fromRoots<C = number>(
    roots: readonly unknown[],
    field: PolynomialField<C> = activeField as PolynomialField<C>,
  ): Polynomial<C> {
    if (!Array.isArray(roots)) {
      throw new TypeError('Polynomial.fromRoots expects an array');
    }

    const zero = field.parse(0);
    const nonZeroRoots = roots.filter(root => !field.equals(field.parse(root), zero));
    const zeroCount = roots.length - nonZeroRoots.length;
    const one = new Polynomial<C>(1, field);

    const product = (values: readonly unknown[]): Polynomial<C> => {
      if (values.length === 0) {
        return one;
      }
      if (values.length === 1) {
        return new Polynomial<C>({
          0: field.mul(field.parse(values[0]), -1),
          1: field.parse(1),
        }, field);
      }

      const middle = Math.floor(values.length / 2);
      return product(values.slice(0, middle)).mul(product(values.slice(middle)));
    };

    const base = product(nonZeroRoots);
    const coefficients = createCoefficientMap<C>();
    for (const [key, coefficient] of Object.entries(base.coeff)) {
      coefficients[Number(key) + zeroCount] = coefficient;
    }
    return new Polynomial<C>(coefficients, field);
  }

  static quotient(
    modulus: PolynomialInput<number>,
  ): PolynomialQuotientRing<number>;
  static quotient<C>(
    modulus: Polynomial<C>,
  ): PolynomialQuotientRing<C>;
  static quotient<C>(
    modulus: PolynomialInput<C>,
    field: PolynomialField<C>,
  ): PolynomialQuotientRing<C>;
  static quotient<C = number>(
    modulus: PolynomialInput<C>,
    field?: PolynomialField<C>,
  ): PolynomialQuotientRing<C> {
    const coefficientField = modulus instanceof Polynomial && field === undefined
      ? modulus.field
      : field ?? activeField as PolynomialField<C>;
    const polynomial = modulus instanceof Polynomial && field === undefined
      ? modulus as Polynomial<C>
      : new Polynomial<C>(
        modulus,
        coefficientField,
      );
    return new QuotientRing(polynomial, coefficientField);
  }

  static enumerate(maxDegree: number, maxElements?: number): Polynomial<number>[];
  static enumerate<C>(
    maxDegree: number,
    field: PolynomialField<C>,
    maxElements?: number,
  ): Polynomial<C>[];
  static enumerate<C = number>(
    maxDegree: number,
    fieldOrLimit?: PolynomialField<C> | number,
    inputLimit?: number,
  ): Polynomial<C>[] {
    if (!Number.isSafeInteger(maxDegree) || maxDegree < 0) {
      throw new RangeError('Maximum polynomial degree must be a non-negative safe integer');
    }

    const field = typeof fieldOrLimit === 'object'
      ? fieldOrLimit
      : activeField as PolynomialField<C>;
    const maxElements = typeof fieldOrLimit === 'number'
      ? fieldOrLimit
      : inputLimit ?? 1_000_000;
    if (!Number.isSafeInteger(maxElements) || maxElements < 0) {
      throw new RangeError('Enumeration limit must be a non-negative safe integer');
    }

    let polynomial = new Polynomial<C>(0, field);
    const characteristic = polynomial.assertPrimeField();
    const count = characteristic ** (maxDegree + 1);
    if (!Number.isSafeInteger(count)) {
      throw new RangeError('Polynomial enumeration size exceeds the safe integer range');
    }
    if (count > maxElements) {
      throw new RangeError(`Polynomial enumeration has ${count} elements, exceeding the limit ${maxElements}`);
    }

    const polynomials: Polynomial<C>[] = [];
    for (let index = 0; index < count; index += 1) {
      polynomials.push(polynomial);
      polynomial = polynomial.successor();
    }
    return polynomials;
  }

  static cyclicCodeGenerators(
    length: number,
    dimension: number,
    maxGenerators?: number,
  ): Polynomial<number>[];
  static cyclicCodeGenerators<C>(
    length: number,
    dimension: number,
    field: PolynomialField<C>,
    maxGenerators?: number,
  ): Polynomial<C>[];
  static cyclicCodeGenerators<C = number>(
    length: number,
    dimension: number,
    fieldOrLimit?: PolynomialField<C> | number,
    inputLimit?: number,
  ): Polynomial<C>[] {
    if (!Number.isSafeInteger(length) || length <= 0) {
      throw new RangeError('Cyclic code length must be a positive safe integer');
    }
    if (!Number.isSafeInteger(dimension) || dimension < 0 || dimension > length) {
      throw new RangeError('Cyclic code dimension must be between zero and the code length');
    }

    const field = typeof fieldOrLimit === 'object'
      ? fieldOrLimit
      : activeField as PolynomialField<C>;
    const maxGenerators = typeof fieldOrLimit === 'number'
      ? fieldOrLimit
      : inputLimit ?? 1_000_000;
    if (!Number.isSafeInteger(maxGenerators) || maxGenerators < 0) {
      throw new RangeError('Generator limit must be a non-negative safe integer');
    }

    const one = field.parse(1);
    const ambient = new Polynomial<C>({
      0: field.mul(one, -1),
      [length]: one,
    }, field);
    const groupedFactors: Array<{ factor: Polynomial<C>; multiplicity: number }> = [];
    for (const factor of ambient.factor()) {
      const group = groupedFactors.find(candidate => candidate.factor.equals(factor));
      if (group === undefined) {
        groupedFactors.push({ factor, multiplicity: 1 });
      } else {
        group.multiplicity += 1;
      }
    }

    const generators: Polynomial<C>[] = [];
    const targetDegree = length - dimension;
    const collect = (
      factorIndex: number,
      remainingDegree: number,
      product: Polynomial<C>,
    ): void => {
      if (factorIndex === groupedFactors.length) {
        if (remainingDegree === 0) {
          if (generators.length >= maxGenerators) {
            throw new RangeError(`Cyclic code generator count exceeds the limit ${maxGenerators}`);
          }
          generators.push(product);
        }
        return;
      }

      const group = groupedFactors[factorIndex]!;
      let factorPower = new Polynomial<C>(1, field);
      for (
        let exponent = 0;
        exponent <= group.multiplicity && exponent * group.factor.degree() <= remainingDegree;
        exponent += 1
      ) {
        collect(
          factorIndex + 1,
          remainingDegree - exponent * group.factor.degree(),
          product.mul(factorPower),
        );
        factorPower = factorPower.mul(group.factor);
      }
    };

    collect(0, targetDegree, new Polynomial<C>(1, field));
    return generators.sort((left, right) => left.toString().localeCompare(right.toString()));
  }

  static setField(field: BuiltInField | PolynomialField<unknown>): void {
    if (typeof field === 'object' && field !== null) {
      assertField(field);
      activeField = field;
      return;
    }
    if (field === 'R') {
      activeField = realField as PolynomialField<unknown>;
      return;
    }

    const constructors: Partial<Record<BuiltInField, AlgebraicConstructor>> = {
      Q: Fraction as unknown as AlgebraicConstructor,
      C: Complex as unknown as AlgebraicConstructor,
      H: Quaternion as unknown as AlgebraicConstructor,
    };
    const constructor = constructors[field];
    if (constructor !== undefined) {
      activeField = createAlgebraicField(constructor) as PolynomialField<unknown>;
      return;
    }

    if (typeof field === 'string' && /^Z\d+$/.test(field)) {
      activeField = createModularField(Number(field.slice(1))) as PolynomialField<unknown>;
      return;
    }
    throw new RangeError(`Unsupported polynomial field: ${String(field)}`);
  }

  private create(coefficients: CoefficientMap<C>): Polynomial<C> {
    return new Polynomial(coefficients, this.field);
  }

  private assertPrimeField(): number {
    const characteristic = this.field.characteristic;
    if (characteristic === undefined || !isPrime(characteristic)) {
      throw new TypeError('Polynomial factorization requires a prime coefficient field');
    }
    return characteristic;
  }

  private coefficientDigit(exponent: number): number {
    const characteristic = this.assertPrimeField();
    const coefficient = this.coeff[exponent] ?? this.field.parse(0);
    for (let digit = 0; digit < characteristic; digit += 1) {
      if (this.field.equals(coefficient, digit)) {
        return digit;
      }
    }
    throw new TypeError('Coefficient is not represented by the prime subfield');
  }

  private factorMonicFinite(polynomial: Polynomial<C>): Polynomial<C>[] {
    if (polynomial.degree() <= 1) {
      return [polynomial.monic()];
    }

    const characteristic = this.assertPrimeField();
    const derivative = polynomial.derive();
    if (derivative.isZero()) {
      const rootCoefficients = createCoefficientMap<C>();
      for (const [key, coefficient] of Object.entries(polynomial.coeff)) {
        const exponent = Number(key);
        if (exponent % characteristic !== 0) {
          throw new Error('Unable to extract the Frobenius root');
        }
        rootCoefficients[exponent / characteristic] = coefficient;
      }

      const rootFactors = this.factorMonicFinite(this.create(rootCoefficients));
      return rootFactors.flatMap(factor => (
        Array.from({ length: characteristic }, () => factor.clone())
      ));
    }

    const repeatedPart = polynomial.gcd(derivative);
    const squareFreePart = polynomial.div(repeatedPart).monic();
    const factors = this.factorSquareFree(squareFreePart);
    if (repeatedPart.degree() > 0) {
      factors.push(...this.factorMonicFinite(repeatedPart.monic()));
    }
    return factors;
  }

  private factorSquareFree(polynomial: Polynomial<C>): Polynomial<C>[] {
    if (polynomial.degree() <= 1) {
      return [polynomial.monic()];
    }

    const characteristic = this.assertPrimeField();
    const basis = this.berlekampBasis(polynomial);
    if (basis.length === 1) {
      return [polynomial.monic()];
    }

    for (const vector of basis) {
      const splitter = this.create(
        Object.fromEntries(vector.map((coefficient, exponent) => [exponent, coefficient])),
      );
      if (splitter.degree() <= 0) {
        continue;
      }

      for (let scalar = 0; scalar < characteristic; scalar += 1) {
        const divisor = polynomial.gcd(splitter.sub(scalar));
        if (divisor.degree() > 0 && divisor.degree() < polynomial.degree()) {
          return [
            ...this.factorSquareFree(divisor.monic()),
            ...this.factorSquareFree(polynomial.div(divisor).monic()),
          ];
        }
      }
    }

    throw new Error('Berlekamp factorization failed to split a reducible polynomial');
  }

  private berlekampBasis(polynomial: Polynomial<C>): C[][] {
    const characteristic = this.assertPrimeField();
    const degree = polynomial.degree();
    const ring = Polynomial.quotient(polynomial);
    const frobenius = ring.create({ 1: this.field.parse(1) }).pow(characteristic);
    let power = ring.one();
    const matrix = Array.from(
      { length: degree },
      () => Array.from({ length: degree }, () => this.field.parse(0)),
    );

    for (let column = 0; column < degree; column += 1) {
      const coefficients = power.value.toArray(degree);
      for (let row = 0; row < degree; row += 1) {
        matrix[row]![column] = this.field.sub(
          coefficients[row] as C,
          row === column ? 1 : 0,
        );
      }
      power = power.mul(frobenius);
    }

    return this.nullspace(matrix);
  }

  private nullspace(matrix: C[][]): C[][] {
    const rowCount = matrix.length;
    const columnCount = matrix[0]?.length ?? 0;
    const reduced = matrix.map(row => [...row]);
    const pivotColumns: number[] = [];
    let pivotRow = 0;

    for (let column = 0; column < columnCount && pivotRow < rowCount; column += 1) {
      let selectedRow = pivotRow;
      while (
        selectedRow < rowCount
        && this.field.empty(reduced[selectedRow]?.[column])
      ) {
        selectedRow += 1;
      }
      if (selectedRow === rowCount) {
        continue;
      }

      [reduced[pivotRow], reduced[selectedRow]] = [
        reduced[selectedRow] as C[],
        reduced[pivotRow] as C[],
      ];
      const pivot = reduced[pivotRow]?.[column] as C;
      for (let index = 0; index < columnCount; index += 1) {
        reduced[pivotRow]![index] = this.field.div(reduced[pivotRow]![index] as C, pivot);
      }

      for (let row = 0; row < rowCount; row += 1) {
        if (row === pivotRow || this.field.empty(reduced[row]?.[column])) {
          continue;
        }
        const factor = reduced[row]?.[column] as C;
        for (let index = 0; index < columnCount; index += 1) {
          reduced[row]![index] = this.field.sub(
            reduced[row]![index] as C,
            this.field.mul(factor, reduced[pivotRow]![index]),
          );
        }
      }

      pivotColumns.push(column);
      pivotRow += 1;
    }

    const pivots = new Set(pivotColumns);
    const freeColumns = Array.from(
      { length: columnCount },
      (_, column) => column,
    ).filter(column => !pivots.has(column));

    return freeColumns.map(freeColumn => {
      const vector = Array.from(
        { length: columnCount },
        () => this.field.parse(0),
      );
      vector[freeColumn] = this.field.parse(1);
      for (let row = 0; row < pivotColumns.length; row += 1) {
        vector[pivotColumns[row] as number] = this.field.mul(
          reduced[row]?.[freeColumn] as C,
          -1,
        );
      }
      return vector;
    });
  }

  private shift(offset: number): Polynomial<C> {
    const result = createCoefficientMap<C>();
    const parsedOffset = this.field.parse(offset);

    for (const [key, coefficient] of Object.entries(this.coeff)) {
      const exponent = Number(key);
      for (let targetExponent = 0; targetExponent <= exponent; targetExponent += 1) {
        const scalar = binomialCoefficient(exponent, targetExponent);
        const offsetPower = this.field.pow(parsedOffset, exponent - targetExponent);
        const term = this.field.mul(this.field.mul(coefficient, scalar), offsetPower);
        const current = result[targetExponent];
        result[targetExponent] = current === undefined
          ? term
          : this.field.add(current, term);
      }
    }

    return this.create(normalize(result, this.field));
  }

  private isSquareRoot(candidate: Polynomial<C>): boolean {
    const square = candidate.mul(candidate);
    const degree = Math.max(this.degree(), square.degree());
    const zero = this.field.parse(0);
    for (let exponent = 0; exponent <= degree; exponent += 1) {
      const expected = this.coeff[exponent] ?? zero;
      const actual = square.coeff[exponent] ?? zero;
      if (!this.field.equals(actual, expected)) {
        if (this.field.characteristic !== 0 || this.field.magnitude === undefined) {
          return false;
        }

        const difference = this.field.magnitude(this.field.sub(actual, expected));
        const scale = Math.max(
          1,
          this.field.magnitude(actual),
          this.field.magnitude(expected),
        );
        if (difference > Number.EPSILON * 32 * scale) {
          return false;
        }
      }
    }
    return true;
  }

  private divideCoefficients(
    numerator: MutableCoefficientMap<C>,
    denominator: CoefficientMap<C>,
  ): MutableCoefficientMap<C> {
    const quotient = createCoefficientMap<C>();
    const denominatorDegree = degreeOf(denominator, this.field);
    if (denominatorDegree < 0) {
      throw new RangeError('Division by the zero polynomial');
    }

    let numeratorDegree = degreeOf(numerator, this.field);
    const trace: Polynomial<C>[] = [];
    while (numeratorDegree >= denominatorDegree) {
      const exponent = numeratorDegree - denominatorDegree;
      const coefficient = this.field.div(
        numerator[numeratorDegree] as C,
        denominator[denominatorDegree] as C,
      );
      quotient[exponent] = coefficient;

      const traceCoefficients = createCoefficientMap<C>();
      for (const [key, denominatorCoefficient] of Object.entries(denominator)) {
        const targetExponent = Number(key) + exponent;
        const product = this.field.mul(denominatorCoefficient, coefficient);
        const current = numerator[targetExponent];
        numerator[targetExponent] = current === undefined
          ? this.field.mul(product, -1)
          : this.field.sub(current, product);
        traceCoefficients[targetExponent] = product;
      }
      normalize(numerator, this.field);

      if (Polynomial.trace !== null) {
        trace.push(this.create(traceCoefficients));
      }
      const nextDegree = degreeOf(numerator, this.field);
      if (nextDegree >= numeratorDegree) {
        throw new Error('Polynomial division did not reduce the remainder degree');
      }
      numeratorDegree = nextDegree;
    }

    if (Polynomial.trace !== null) {
      trace.push(this.create(numerator));
      Polynomial.trace = trace as Polynomial<unknown>[];
    }
    return normalize(quotient, this.field);
  }

  private format(method: 'toString' | 'toLatex'): string {
    const exponents = Object.keys(this.coeff).map(Number).sort((left, right) => right - left);
    let result = '';

    for (const exponent of exponents) {
      const coefficient = this.coeff[exponent];
      if (coefficient === undefined || coefficient === null) {
        continue;
      }

      if (isComplex(coefficient)) {
        result += this.formatComplexTerm(coefficient, exponent, result.length > 0);
        continue;
      }

      const value = (coefficient as { valueOf?: () => unknown }).valueOf?.() ?? coefficient;
      if (value === 0) {
        continue;
      }
      if (result.length > 0 && Number(value) > 0) {
        result += '+';
      }
      if (value === -1 && exponent !== 0) {
        result += '-';
      } else if (value !== 1 || exponent === 0) {
        const formatter = (coefficient as unknown as Record<string, unknown>)[method];
        result += typeof formatter === 'function'
          ? String(formatter.call(coefficient))
          : String(coefficient);
      }
      result += this.formatVariable(exponent);
    }
    return result || '0';
  }

  private formatComplexTerm(value: Complex, exponent: number, hasPrevious: boolean): string {
    let result = '';
    if (value.re !== 0) {
      if (hasPrevious && value.re > 0) {
        result += '+';
      }
      if (value.re === -1 && exponent !== 0) {
        result += '-';
      } else if (value.re !== 1 || exponent === 0) {
        result += String(value.re);
      }
      result += this.formatVariable(exponent);
      hasPrevious = true;
    }

    if (value.im !== 0) {
      if (hasPrevious && value.im > 0) {
        result += '+';
      }
      if (value.im === -1) {
        result += '-';
      } else if (value.im !== 1) {
        result += String(value.im);
      }
      result += `i${this.formatVariable(exponent)}`;
    }
    return result;
  }

  private formatVariable(exponent: number): string {
    if (exponent === 0) {
      return '';
    }
    if (exponent === 1) {
      return 'x';
    }
    return `x^${exponent}`;
  }
}

class QuotientRing<C> implements PolynomialQuotientRing<C> {
  readonly modulus: Polynomial<C>;
  readonly field: PolynomialField<C>;
  private fieldStatus?: boolean;

  constructor(modulus: Polynomial<C>, field: PolynomialField<C>) {
    if (modulus.degree() < 1) {
      throw new RangeError('Quotient ring modulus must have positive degree');
    }
    this.field = field;
    this.modulus = modulus.monic();
  }

  create(input?: PolynomialResidueInput<C>): PolynomialResidue<C> {
    if (isPolynomialResidue(input)) {
      const residue = input as Residue<C>;
      this.assertCompatible(residue);
      return new Residue(this, residue.value);
    }
    return new Residue(this, this.polynomial(input));
  }

  zero(): PolynomialResidue<C> {
    return this.create(0);
  }

  one(): PolynomialResidue<C> {
    return this.create(1);
  }

  evaluate(
    input: PolynomialInput<C>,
    value: PolynomialResidueInput<C>,
  ): PolynomialResidue<C> {
    const polynomial = this.polynomial(input);
    const element = this.create(value);
    let result = this.zero();

    for (let exponent = polynomial.degree(); exponent >= 0; exponent -= 1) {
      result = result.mul(element);
      const coefficient = polynomial.coeff[exponent];
      if (!this.field.empty(coefficient)) {
        result = result.add({ 0: coefficient as C });
      }
    }
    return result;
  }

  isRoot(
    polynomial: PolynomialInput<C>,
    value: PolynomialResidueInput<C>,
  ): boolean {
    return this.evaluate(polynomial, value).isZero();
  }

  order(): number {
    const characteristic = this.assertPrimeCharacteristic();
    const order = characteristic ** this.modulus.degree();
    if (!Number.isSafeInteger(order)) {
      throw new RangeError('Quotient ring order exceeds the safe integer range');
    }
    return order;
  }

  isField(): boolean {
    this.assertPrimeCharacteristic();
    this.fieldStatus ??= this.modulus.isIrreducible();
    return this.fieldStatus;
  }

  elements(maxElements = 1_000_000): PolynomialResidue<C>[] {
    const order = this.order();
    this.assertEnumerationLimit(maxElements, order);
    const characteristic = this.assertPrimeCharacteristic();
    const degree = this.modulus.degree();
    const elements: PolynomialResidue<C>[] = [];

    for (let index = 0; index < order; index += 1) {
      let encoded = index;
      const coefficients = createCoefficientMap<C>();
      for (let exponent = 0; exponent < degree; exponent += 1) {
        const digit = encoded % characteristic;
        if (digit !== 0) {
          coefficients[exponent] = this.field.parse(digit);
        }
        encoded = Math.floor(encoded / characteristic);
      }
      elements.push(this.create(coefficients));
    }
    return elements;
  }

  primitiveElements(maxElements = 1_000_000): PolynomialResidue<C>[] {
    this.assertFiniteField();
    return this.elements(maxElements).filter(element => element.isPrimitive());
  }

  primitivePolynomials(maxElements = 1_000_000): Polynomial<C>[] {
    const polynomials: Polynomial<C>[] = [];
    for (const element of this.primitiveElements(maxElements)) {
      const minimalPolynomial = element.minimalPolynomial();
      if (!polynomials.some(polynomial => polynomial.equals(minimalPolynomial))) {
        polynomials.push(minimalPolynomial);
      }
    }
    return polynomials.sort((left, right) => left.compare(right));
  }

  polynomial(input?: PolynomialInput<C>): Polynomial<C> {
    return this.modulus.mul(0).add(input);
  }

  assertCompatible(input: Residue<C>): void {
    if (input.ring !== this) {
      throw new RangeError('Cannot combine residue classes from different quotient rings');
    }
  }

  assertFiniteField(): void {
    if (!this.isField()) {
      throw new RangeError('Operation requires an irreducible quotient-ring modulus');
    }
  }

  private assertPrimeCharacteristic(): number {
    const characteristic = this.field.characteristic;
    if (characteristic === undefined || !isPrime(characteristic)) {
      throw new TypeError('Operation requires a finite prime coefficient field');
    }
    return characteristic;
  }

  private assertEnumerationLimit(maxElements: number, order: number): void {
    if (!Number.isSafeInteger(maxElements) || maxElements < 0) {
      throw new RangeError('Enumeration limit must be a non-negative safe integer');
    }
    if (order > maxElements) {
      throw new RangeError(`Quotient ring has ${order} elements, exceeding the limit ${maxElements}`);
    }
  }
}

class Residue<C> implements PolynomialResidue<C> {
  readonly ring: QuotientRing<C>;
  readonly value: Polynomial<C>;

  constructor(ring: QuotientRing<C>, value: Polynomial<C>) {
    this.ring = ring;
    this.value = value.mod(ring.modulus);
  }

  add(input: PolynomialResidueInput<C>): PolynomialResidue<C> {
    return this.ring.create(this.value.add(this.valueOf(input)));
  }

  sub(input: PolynomialResidueInput<C>): PolynomialResidue<C> {
    return this.ring.create(this.value.sub(this.valueOf(input)));
  }

  mul(input: PolynomialResidueInput<C>): PolynomialResidue<C> {
    return this.ring.create(this.value.mul(this.valueOf(input)));
  }

  div(input: PolynomialResidueInput<C>): PolynomialResidue<C> {
    const divisor = this.ring.create(input);
    return this.mul(divisor.inverse());
  }

  neg(): PolynomialResidue<C> {
    return this.ring.create(this.value.neg());
  }

  pow(exponent: number): PolynomialResidue<C> {
    if (!Number.isSafeInteger(exponent) || exponent < 0) {
      throw new RangeError('Residue exponent must be a non-negative safe integer');
    }

    let result = this.ring.one();
    let factor: PolynomialResidue<C> = this;
    let power = exponent;
    while (power > 0) {
      if (power % 2 === 1) {
        result = result.mul(factor);
      }
      power = Math.floor(power / 2);
      if (power > 0) {
        factor = factor.mul(factor);
      }
    }
    return result;
  }

  inverse(): PolynomialResidue<C> {
    if (this.isZero()) {
      throw new RangeError('Zero residue class is not invertible');
    }

    let oldRemainder = this.ring.modulus;
    let remainder = this.value;
    let oldCoefficient = this.ring.polynomial(0);
    let coefficient = this.ring.polynomial(1);

    while (!remainder.isZero()) {
      const quotient = oldRemainder.div(remainder);
      [oldRemainder, remainder] = [
        remainder,
        oldRemainder.sub(quotient.mul(remainder)),
      ];
      [oldCoefficient, coefficient] = [
        coefficient,
        oldCoefficient.sub(quotient.mul(coefficient)),
      ];
    }

    if (oldRemainder.degree() !== 0) {
      throw new RangeError('Residue class is not invertible');
    }
    return this.ring.create(oldCoefficient.div(oldRemainder.lc() as C));
  }

  additiveOrder(): number {
    const characteristic = this.ring.field.characteristic;
    if (characteristic === undefined || !isPrime(characteristic)) {
      throw new TypeError('Additive order requires a finite prime coefficient field');
    }
    return this.isZero() ? 1 : characteristic;
  }

  additiveSubgroup(maxElements = 1_000_000): PolynomialResidue<C>[] {
    const order = this.additiveOrder();
    this.assertSubgroupLimit(maxElements, order);
    const subgroup: PolynomialResidue<C>[] = [this.ring.zero()];
    let element = this.ring.zero();
    for (let index = 1; index < order; index += 1) {
      element = element.add(this);
      subgroup.push(element);
    }
    return subgroup;
  }

  multiplicativeOrder(): number {
    this.ring.assertFiniteField();
    if (this.isZero()) {
      throw new RangeError('Zero has no multiplicative order');
    }

    let order = this.ring.order() - 1;
    for (const divisor of primeDivisors(order)) {
      while (order % divisor === 0 && this.pow(order / divisor).equals(1)) {
        order /= divisor;
      }
    }
    return order;
  }

  multiplicativeSubgroup(maxElements = 1_000_000): PolynomialResidue<C>[] {
    const order = this.multiplicativeOrder();
    this.assertSubgroupLimit(maxElements, order);
    const subgroup: PolynomialResidue<C>[] = [this.ring.one()];
    let element = this.ring.one();
    for (let exponent = 1; exponent < order; exponent += 1) {
      element = element.mul(this);
      subgroup.push(element);
    }
    return subgroup;
  }

  isPrimitive(): boolean {
    this.ring.assertFiniteField();
    return !this.isZero() && this.multiplicativeOrder() === this.ring.order() - 1;
  }

  minimalPolynomial(): Polynomial<C> {
    this.ring.assertFiniteField();
    const characteristic = this.ring.field.characteristic as number;
    const conjugates: PolynomialResidue<C>[] = [];
    let conjugate: PolynomialResidue<C> = this;

    do {
      conjugates.push(conjugate);
      conjugate = conjugate.pow(characteristic);
      if (conjugates.length > this.ring.modulus.degree()) {
        throw new Error('Frobenius orbit exceeds the extension degree');
      }
    } while (!conjugate.equals(this));

    let coefficients: PolynomialResidue<C>[] = [this.ring.one()];
    for (const root of conjugates) {
      const next = Array.from(
        { length: coefficients.length + 1 },
        () => this.ring.zero(),
      );
      for (let exponent = 0; exponent < coefficients.length; exponent += 1) {
        const coefficient = coefficients[exponent] as PolynomialResidue<C>;
        next[exponent] = next[exponent]!.sub(coefficient.mul(root));
        next[exponent + 1] = next[exponent + 1]!.add(coefficient);
      }
      coefficients = next;
    }

    const result = createCoefficientMap<C>();
    for (let exponent = 0; exponent < coefficients.length; exponent += 1) {
      const representative = coefficients[exponent]!.value;
      if (representative.degree() > 0) {
        throw new Error('Minimal polynomial coefficient is not in the base field');
      }
      const coefficient = representative.coeff[0] ?? this.ring.field.parse(0);
      if (!this.ring.field.empty(coefficient)) {
        result[exponent] = coefficient;
      }
    }
    return this.ring.polynomial(result);
  }

  equals(input: PolynomialResidueInput<C>): boolean {
    return this.value.equals(this.valueOf(input));
  }

  isZero(): boolean {
    return this.value.isZero();
  }

  private assertSubgroupLimit(maxElements: number, order: number): void {
    if (!Number.isSafeInteger(maxElements) || maxElements < 0) {
      throw new RangeError('Subgroup limit must be a non-negative safe integer');
    }
    if (order > maxElements) {
      throw new RangeError(`Subgroup has ${order} elements, exceeding the limit ${maxElements}`);
    }
  }

  toString(): string {
    return this.value.toString();
  }

  private valueOf(input: PolynomialResidueInput<C>): Polynomial<C> {
    if (isPolynomialResidue(input)) {
      const residue = input as Residue<C>;
      this.ring.assertCompatible(residue);
      return residue.value;
    }
    return this.ring.polynomial(input);
  }
}

function isPolynomialResidue<C>(
  input: PolynomialResidueInput<C>,
): input is PolynomialResidue<C> {
  return input instanceof Residue;
}

export default Polynomial;