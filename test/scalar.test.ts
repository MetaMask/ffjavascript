import { describe, expect, it } from "vitest";

import * as Scalar from "../src/scalar";

describe("Basic scalar conversions", () => {
  it("Should convert Native values", () => {
    expect(Scalar.eq(Scalar.e("0x12")!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e("0x12", 16)!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e("12", 16)!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e("18")!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e("18", 10)!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e(18, 10)!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e(18n, 10)!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e(0x12, 10)!, 18)).toBe(true);
    expect(Scalar.eq(Scalar.e(0x12n, 10)!, 18)).toBe(true);
  });

  it("Should convert to js Number Native", () => {
    const maxJsNum = Number.MAX_SAFE_INTEGER;
    const maxToScalar = Scalar.e(maxJsNum)!;

    const backToNum = Scalar.toNumber(maxToScalar);
    expect(backToNum).toBe(maxJsNum);

    const overMaxJsNum = Scalar.add(maxToScalar, 1);
    expect(() => Scalar.toNumber(overMaxJsNum)).toThrow();
  });
});
