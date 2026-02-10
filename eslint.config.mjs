import toruslabsJavascript from "@toruslabs/eslint-config-javascript";

export default [
  ...toruslabsJavascript,
  {
    ignores: ["dist/"],
    rules: {
      "no-redeclare": ["error", { builtinGlobals: false }],
    },
  },
];
