/**
 * @license Polynomial.js v2.0.0 9/10/2026
 * https://github.com/rawify/Polynomial.js
 *
 * Copyright (c) 2025-2026, Robert Eisele (https://raw.org/)
 * Licensed under the MIT license.
 */
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
export type PolynomialInput<C> = Polynomial<C> | CoefficientMap<C> | readonly unknown[] | C | string | number | null | undefined;
export type PolynomialResidueInput<C> = PolynomialResidue<C> | PolynomialInput<C>;
export interface PolynomialQuotientRing<C> {
    readonly modulus: Polynomial<C>;
    create(input?: PolynomialResidueInput<C>): PolynomialResidue<C>;
    zero(): PolynomialResidue<C>;
    one(): PolynomialResidue<C>;
    evaluate(polynomial: PolynomialInput<C>, value: PolynomialResidueInput<C>): PolynomialResidue<C>;
    isRoot(polynomial: PolynomialInput<C>, value: PolynomialResidueInput<C>): boolean;
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
export declare class Polynomial<C = number> {
    static trace: boolean | Polynomial<unknown>[] | null;
    readonly coeff: MutableCoefficientMap<C>;
    private readonly field;
    constructor(input?: PolynomialInput<C>, field?: PolynomialField<C>);
    gcd(input: PolynomialInput<C>): Polynomial<C>;
    lcm(input: PolynomialInput<C>): Polynomial<C>;
    neg(): Polynomial<C>;
    reciprocal(): Polynomial<C>;
    sqrt(): Polynomial<C> | null;
    eval(value: unknown): C | number;
    compose(input: PolynomialInput<C>): Polynomial<C>;
    /** @deprecated Use eval(). */
    result(value: unknown): C | number;
    lc(): C | undefined;
    lm(): Polynomial<C>;
    monic(): Polynomial<C>;
    equals(input: PolynomialInput<C>): boolean;
    isZero(): boolean;
    isMonic(): boolean;
    compare(input: PolynomialInput<C>): -1 | 0 | 1;
    successor(): Polynomial<C>;
    isDivisibleBy(input: PolynomialInput<C>): boolean;
    factor(): Polynomial<C>[];
    isIrreducible(): boolean;
    cyclicCodeBasis(length: number): Polynomial<C>[];
    dispersion(maxShift?: number): number;
    add(input: PolynomialInput<C>): Polynomial<C>;
    sub(input: PolynomialInput<C>): Polynomial<C>;
    mul(input: PolynomialInput<C>): Polynomial<C>;
    addmul(left: PolynomialInput<C>, right: PolynomialInput<C>): Polynomial<C>;
    div(input: PolynomialInput<C>): Polynomial<C>;
    mod(input: PolynomialInput<C>): Polynomial<C>;
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
    static fromRoots<C = number>(roots: readonly unknown[], field?: PolynomialField<C>): Polynomial<C>;
    static quotient(modulus: PolynomialInput<number>): PolynomialQuotientRing<number>;
    static quotient<C>(modulus: Polynomial<C>): PolynomialQuotientRing<C>;
    static quotient<C>(modulus: PolynomialInput<C>, field: PolynomialField<C>): PolynomialQuotientRing<C>;
    static enumerate(maxDegree: number, maxElements?: number): Polynomial<number>[];
    static enumerate<C>(maxDegree: number, field: PolynomialField<C>, maxElements?: number): Polynomial<C>[];
    static cyclicCodeGenerators(length: number, dimension: number, maxGenerators?: number): Polynomial<number>[];
    static cyclicCodeGenerators<C>(length: number, dimension: number, field: PolynomialField<C>, maxGenerators?: number): Polynomial<C>[];
    static setField(field: BuiltInField | PolynomialField<unknown>): void;
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
export default Polynomial;
