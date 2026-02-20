import toruslabsTypescript from "@toruslabs/eslint-config-typescript";

export default [
  ...toruslabsTypescript,
  {
    ignores: ["dist/"],
    rules: {
      "no-redeclare": ["error", { builtinGlobals: false }],
    },
  },
];
