import { defineConfig } from "vitest/config";
import { playwright } from "@vitest/browser-playwright";

export default defineConfig({
  test: {
    reporters: "verbose",
    browser: {
      screenshotFailures: false,
      headless: true,
      provider: playwright(),
      instances: [
        {
          name: "Chrome",
          browser: "chromium",
        },
        { name: "Firefox", browser: "firefox" },
        { name: "Safari", browser: "webkit" },
      ],
      enabled: true,
    },
    coverage: {
      reporter: ["text"],
      provider: "istanbul",
      include: ["src/**/*.js"],
    },
  },
});
