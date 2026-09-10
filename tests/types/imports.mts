import Polynomial, {
  Polynomial as NamedPolynomial,
  type BuiltInField,
  type CoefficientMap,
  type PolynomialField,
  type PolynomialInput,
  type PolynomialQuotientRing,
  type PolynomialResidue,
} from 'polynomial';

const fieldName: BuiltInField = 'Z7';
const coefficients: CoefficientMap<number> = { 0: 1, 2: 3 };
const input: PolynomialInput<number> = coefficients;
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

const polynomial: Polynomial<number> = new NamedPolynomial(input, field);
const fromRoots: Polynomial<number> = Polynomial.fromRoots([1, 2], field);
const dispersion: number = polynomial.dispersion();
const squareRoot: Polynomial<number> | null = polynomial.sqrt();
const composition: Polynomial<number> = polynomial.compose('x+1');
const leastCommonMultiple: Polynomial<number> = polynomial.lcm('x-1');
const coefficientsArray: number[] = polynomial.toArray(4);
const equal: boolean = polynomial.equals('3x^2+1');
const zero: boolean = polynomial.isZero();
const monic: boolean = polynomial.isMonic();
const divisible: boolean = polynomial.isDivisibleBy('x+1');
const quotientRing: PolynomialQuotientRing<number> = Polynomial.quotient('x^2+1', field);
const residue: PolynomialResidue<number> = quotientRing.create('x');
const evaluatedResidue: PolynomialResidue<number> = quotientRing.evaluate(polynomial, residue);
const residueRoot: boolean = quotientRing.isRoot(polynomial, residue);
const inverseResidue: PolynomialResidue<number> = residue.inverse();
const factors: Polynomial<number>[] = polynomial.factor();
const irreducible: boolean = polynomial.isIrreducible();
const fieldOrder: number = quotientRing.order();
const fieldElements: PolynomialResidue<number>[] = quotientRing.elements();
const primitiveElements: PolynomialResidue<number>[] = quotientRing.primitiveElements();
const primitivePolynomials: Polynomial<number>[] = quotientRing.primitivePolynomials();
const additiveOrder: number = residue.additiveOrder();
const additiveSubgroup: PolynomialResidue<number>[] = residue.additiveSubgroup();
const multiplicativeOrder: number = residue.multiplicativeOrder();
const multiplicativeSubgroup: PolynomialResidue<number>[] = residue.multiplicativeSubgroup();
const primitive: boolean = residue.isPrimitive();
const minimalPolynomial: Polynomial<number> = residue.minimalPolynomial();
const cyclicBasis: Polynomial<number>[] = polynomial.cyclicCodeBasis(3);
const cyclicGenerators: Polynomial<number>[] = Polynomial.cyclicCodeGenerators(3, 1, field);
const comparison: -1 | 0 | 1 = polynomial.compare('x');
const successor: Polynomial<number> = polynomial.successor();
const enumeration: Polynomial<number>[] = Polynomial.enumerate(2, field);

void fieldName;
void polynomial;
void fromRoots;
void dispersion;
void squareRoot;
void composition;
void leastCommonMultiple;
void coefficientsArray;
void equal;
void zero;
void monic;
void divisible;
void quotientRing;
void residue;
void evaluatedResidue;
void residueRoot;
void inverseResidue;
void factors;
void irreducible;
void fieldOrder;
void fieldElements;
void primitiveElements;
void primitivePolynomials;
void additiveOrder;
void additiveSubgroup;
void multiplicativeOrder;
void multiplicativeSubgroup;
void primitive;
void minimalPolynomial;
void cyclicBasis;
void cyclicGenerators;
void comparison;
void successor;
void enumeration;