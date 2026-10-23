# Fiji preparation and CI

Run from a checkout, in the installed UMA environment:

```bash
python scripts/prepare_fiji.py --check
```

The supported toolchain is Python 3.10, OpenJDK 11, Maven 3.9.9, and the
unchanged endpoint `sc.fiji:fiji:2.14.0`. The helper registers one historical
dependency before jgo resolves Fiji. Maven uses its configured local repository
(normally `~/.m2/repository`); use the same environment and Maven configuration
for preparation and image analysis. No administrator access is required.

## Artifact provenance

| Item | Value |
|---|---|
| Maven coordinates | `woolz:JWlz:1.4.0` |
| Primary source | [Official Fiji archive](https://sites.imagej.net/Fiji/jars/JWlz-1.4.0.jar-20140304121539) |
| Fallback | [Official Bio-Formats archive](https://downloads.openmicroscopy.org/bio-formats/5.0.0-rc2/artifacts/JWlz-1.4.0.jar) |
| SHA-256 | `0662cfe9600c24cbef591e17f0bfbbbebbae8b0e8d83098735454a999166410f` |
| Verified size | 117846 bytes |
| Install plugin | `org.apache.maven.plugins:maven-install-plugin:3.1.3` |

On 2026-10-23 the SciJava artifact and POM URLs returned HTTP 404. Both official
archive copies had the same digest. The helper checks that digest before any
installation, then generates a coordinate-only POM using Maven's install
plugin. It does not replace the complete Fiji distribution or exclude readers.

Downloads are cached in `~/.cache/uma-tools/bootstrap`. Override this directory
with `--artifact-cache PATH`. A changed or truncated JAR is rejected, including
in the cache; remove only that cached `JWlz-1.4.0.jar` and rerun to download it
again. Do not change the expected digest merely to silence an error. TLS checks
remain enabled. Temporary HTTP failures and timeouts receive at most three
attempts per source; HTTP 404 is not retried at the same URL.

`--check` runs Fiji and a Gaussian filter in a child process, requires Java 11,
disposes ImageJ, closes its workers, and requires process exit within 15 minutes
per attempt. A reported temporary Maven bootstrap error (HTTP 408/429/5xx or
connection/read timeout/reset) permits at most three fresh-process attempts.
Retries pass `-U` through `MAVEN_ARGS` so Maven rechecks failed downloads.
Missing artifacts, wrong filter output, Java mismatches and process timeouts
fail immediately; scientific test failures are never retried by this helper. Running without `--check` only prepares
the dependency. Include the full Maven output when reporting startup failures.

## Independent CI results

- Environment construction and imports: `build-env.yml`.
- Package installation and lightweight command entry points: `core-assays.yml`.
- Required lint and format checks for main, functional, and nuclei packages:
  `lint.yml`, using the same `requirements-dev.txt` as local development.
- Python-only main and functional tests: separate jobs in `test.yml`.
- Real Fiji, calibration, channel equivalence, partial recovery, and process
  exit: dedicated Linux/macOS jobs in `test.yml`. In-process JVM suites use
  separate Python interpreters because a JVM cannot be restarted in place.

Fiji caches are saved only after a successful native preflight. Keys include
OS, architecture, environment, package requirements, Fiji endpoint and the
preparation helper. No broad fallback cache keys are used. The manual
`cold_cache` input skips cache restoration; the first run for a new key also
starts without a restored cache. Each job repeats preparation and JVM checks
in a new process to verify reuse. Bootstrap logs and JUnit results are uploaded.

The older `alternative_assays` and `original_fibronectin_alignment_analysis`
scripts have an explicitly advisory lint audit. Its actual outcome and full
findings are retained in the job summary and `archive-style-findings` artifact.
These findings have not been fixed or silently suppressed. All active-package
style checks and scientific tests remain mandatory pass/fail jobs.
