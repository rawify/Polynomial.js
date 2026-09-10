/**
 * @license Polynomial.js v2.0.0 9/10/2026
 * https://github.com/rawify/Polynomial.js
 *
 * Copyright (c) 2025-2026, Robert Eisele (https://raw.org/)
 * Licensed under the MIT license.
 */
type BuiltInFieldRoot = 'R' | 'Q' | 'C' | 'H' | `Z${number}`;
interface PolynomialFieldRoot<C> {
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
type CoefficientMapRoot<C> = Readonly<Record<number, C>>;
type PolynomialInputRoot<C> = Polynomial<C> | CoefficientMapRoot<C> | readonly unknown[] | C | string | number | null | undefined;
type PolynomialResidueInputRoot<C> = PolynomialResidueRoot<C> | PolynomialInputRoot<C>;
interface PolynomialQuotientRingRoot<C> {
    readonly modulus: Polynomial<C>;
    create(input?: PolynomialResidueInputRoot<C>): PolynomialResidueRoot<C>;
    zero(): PolynomialResidueRoot<C>;
    one(): PolynomialResidueRoot<C>;
    evaluate(polynomial: PolynomialInputRoot<C>, value: PolynomialResidueInputRoot<C>): PolynomialResidueRoot<C>;
    isRoot(polynomial: PolynomialInputRoot<C>, value: PolynomialResidueInputRoot<C>): boolean;
    order(): number;
    isField(): boolean;
    elements(maxElements?: number): PolynomialResidueRoot<C>[];
    primitiveElements(maxElements?: number): PolynomialResidueRoot<C>[];
    primitivePolynomials(maxElements?: number): Polynomial<C>[];
}
interface PolynomialResidueRoot<C> {
    readonly ring: PolynomialQuotientRingRoot<C>;
    readonly value: Polynomial<C>;
    add(input: PolynomialResidueInputRoot<C>): PolynomialResidueRoot<C>;
    sub(input: PolynomialResidueInputRoot<C>): PolynomialResidueRoot<C>;
    mul(input: PolynomialResidueInputRoot<C>): PolynomialResidueRoot<C>;
    div(input: PolynomialResidueInputRoot<C>): PolynomialResidueRoot<C>;
    neg(): PolynomialResidueRoot<C>;
    pow(exponent: number): PolynomialResidueRoot<C>;
    inverse(): PolynomialResidueRoot<C>;
    additiveOrder(): number;
    additiveSubgroup(maxElements?: number): PolynomialResidueRoot<C>[];
    multiplicativeOrder(): number;
    multiplicativeSubgroup(maxElements?: number): PolynomialResidueRoot<C>[];
    isPrimitive(): boolean;
    minimalPolynomial(): Polynomial<C>;
    equals(input: PolynomialResidueInputRoot<C>): boolean;
    isZero(): boolean;
    toString(): string;
}
type MutableCoefficientMapRoot<C> = Record<number, C>;
declare class Polynomial<C = number> {
    static trace: boolean | Polynomial<unknown>[] | null;
    readonly coeff: MutableCoefficientMapRoot<C>;
    private readonly field;
    constructor(input?: PolynomialInputRoot<C>, field?: PolynomialFieldRoot<C>);
    gcd(input: PolynomialInputRoot<C>): Polynomial<C>;
    lcm(input: PolynomialInputRoot<C>): Polynomial<C>;
    neg(): Polynomial<C>;
    reciprocal(): Polynomial<C>;
    sqrt(): Polynomial<C> | null;
    eval(value: unknown): C | number;
    compose(input: PolynomialInputRoot<C>): Polynomial<C>;
    /** @deprecated Use eval(). */
    result(value: unknown): C | number;
    lc(): C | undefined;
    lm(): Polynomial<C>;
    monic(): Polynomial<C>;
    equals(input: PolynomialInputRoot<C>): boolean;
    isZero(): boolean;
    isMonic(): boolean;
    compare(input: PolynomialInputRoot<C>): -1 | 0 | 1;
    successor(): Polynomial<C>;
    isDivisibleBy(input: PolynomialInputRoot<C>): boolean;
    factor(): Polynomial<C>[];
    isIrreducible(): boolean;
    cyclicCodeBasis(length: number): Polynomial<C>[];
    dispersion(maxShift?: number): number;
    add(input: PolynomialInputRoot<C>): Polynomial<C>;
    sub(input: PolynomialInputRoot<C>): Polynomial<C>;
    mul(input: PolynomialInputRoot<C>): Polynomial<C>;
    addmul(left: PolynomialInputRoot<C>, right: PolynomialInputRoot<C>): Polynomial<C>;
    div(input: PolynomialInputRoot<C>): Polynomial<C>;
    mod(input: PolynomialInputRoot<C>): Polynomial<C>;
    pow(exponent: number): Polynomial<C>;
    derive(): Polynomial<C>;
    derive(order: number): Polynomial<C> | null;
    integrate(): Polynomial<C>;
    integrate(order: number): Polynomial<C> | null;
    degree(): number;
    toArray(width?: number): C[];
    clone(): Polynomial<C>;
    toString(): string;
    toLatex(): string;
    toHorner(): string;
    static fromRoots<C = number>(roots: readonly unknown[], field?: PolynomialFieldRoot<C>): Polynomial<C>;
    static quotient(modulus: PolynomialInputRoot<number>): PolynomialQuotientRingRoot<number>;
    static quotient<C>(modulus: Polynomial<C>): PolynomialQuotientRingRoot<C>;
    static quotient<C>(modulus: PolynomialInputRoot<C>, field: PolynomialFieldRoot<C>): PolynomialQuotientRingRoot<C>;
    static enumerate(maxDegree: number, maxElements?: number): Polynomial<number>[];
    static enumerate<C>(maxDegree: number, field: PolynomialFieldRoot<C>, maxElements?: number): Polynomial<C>[];
    static cyclicCodeGenerators(length: number, dimension: number, maxGenerators?: number): Polynomial<number>[];
    static cyclicCodeGenerators<C>(length: number, dimension: number, field: PolynomialFieldRoot<C>, maxGenerators?: number): Polynomial<C>[];
    static setField(field: BuiltInFieldRoot | PolynomialFieldRoot<unknown>): void;
    private create;
    private assertPrimeField;
    private coefficientDigit;
    private factorMonicFinite;
    private factorSquareFree;
    private berlekampBasis;
    private nullspace;
    private shift;
    private isSquareRoot;
    private divideCoefficients;
    private format;
    private formatComplexTerm;
    private formatVariable;
}
declare namespace Polynomial {
  type FieldName = BuiltInFieldRoot;
  type Field<C> = PolynomialFieldRoot<C>;
  type Coefficients<C> = CoefficientMapRoot<C>;
  type Input<C> = PolynomialInputRoot<C>;
  type QuotientRing<C> = PolynomialQuotientRingRoot<C>;
  type Residue<C> = PolynomialResidueRoot<C>;
  type ResidueInput<C> = PolynomialResidueInputRoot<C>;
}
export = Polynomial;
