# rkalos implementation status

Updated **2026-11-24**. Persephone uses the independent rkalos package.
The R package runs without Persephone or Java. Numerical roadmap and canonical reference evidence belong to
[kalos](../kalos/IMPLEMENTATION_STATUS.md).

## M6 public selected-model residual diagnostics in R (2026-11-24)

Version `0.1.0.9018` returns the residual standard error, Ljung–Box
independence and Doornik–Hansen normality results for TR3 preprocessing
and RSA3 fits. The Ljung–Box degrees of freedom subtract the selected
regular and seasonal AR and MA counts from 24. All **65 public JVM
diagnostic values** matched bit-for-bit across 13 raw models with 21,
22 or 23 degrees of freedom. At 24 residuals, independence is undefined
and the normality probability remains available.

Both serial and two-worker OpenMP passed **1,067/1,067** RSA3,
**990/990** RSA1, **130/130** TR1 and **498/498** RSA0 checks and the
isolated three-template Java-free smoke test. The full serial source
suite for `0.1.0.9016` passed **303,509/303,509** across 63 files.
[Selected-model diagnostic replay](inst/tinytest/test-native-tramoseats-rsa3.R)

## Selected-model TR3 residual diagnostics checkpoint (2026-11-24)

The package and CI now pin Kalos
`a4dc93fce2d300a84544f6122e8f231954f30aec`. Its public monthly
residual diagnostics account for the selected regular and seasonal AR
and MA orders. Two JVM captures cover all 13 previously qualified raw
RSA3 models. Offline verification passed for **199 core files, 112
fixtures and 375 evidence files**.

## M5 expanded monthly TR3 and RSA3 identification (2026-11-24)

Version `0.1.0.9017` handles the seasonal controller restart, stationary
scans with the newly qualified orders and a fourth estimation pass after
the fast regression test removes the fitted mean. The R regression
matrix now covers 13 distinct public monthly models. The native orders,
fitted coefficients, likelihoods, **1,196 dated residuals and 8,340
values in 130 component rows** agree with the pinned JVM bit-for-bit.
A rejected first reference comparison still returns its engine error.

Serial and two-worker OpenMP each passed **984/984** RSA3, **990/990**
RSA1, **130/130** standalone TR1 and **498/498** RSA0 checks. Both
installed packages also completed a Java-free smoke fit for each of the
three raw templates. [RSA3 model and component replay](inst/tinytest/test-native-tramoseats-rsa3.R)

## Monthly TR3 seasonal restarts and fourth-pass refits (2026-11-24)

The new core and CI pin `66444bf9262a4a566b53f36fd97c77762e9243d5`
includes the seasonal controller restart, the measured stationary
outlier-scan designs and fourth-pass refitting after the mean is
removed. The two new public captures expand the RSA3 reference
matrix by seven models. Offline verification passed for **199 core
files, 112 fixtures and 372 evidence files**.

## M5 public original-series forecast errors on RSA1 (2026-11-24)

Version `0.1.0.9016` exposes the dated TR1 preprocessing standard errors
in `native_tramoseats_fit(series, "RSA1")$series_forecast_stdev`.
Twelve level and log models matched all **288** values in the public
24-month JVM projection within `2e-8 * (1 + abs(reference))`.
Their 288 distinct SEATS conditional standard errors remained bit-exact.
Seven fits also matched the first six public errors with a custom
3-month backcast and 6-month forecast. A zero horizon returns an empty
vector and the model reports its own numerical availability status.

Serial and two-worker OpenMP each passed **990/990** RSA1,
**448/448** RSA3 and **498/498** RSA0 checks and the isolated Java-free
smoke check. The complete OpenMP source check at `0.1.0.9014` passed
**302,994/302,994** tests, with two known platform compiler warnings.
[RSA1 forecasts and uncertainty](inst/tinytest/test-native-tramoseats-rsa1.R)

## RSA1 original-series forecast error checkpoint (2026-11-24)

Kalos `5a72228cacc78ede80d460a56f9f317c58ddfea5` supplies the
TR1 forecast projection on the original scale, including fitted-regressor
uncertainty and the log-to-level correction. Both R platform builds
compile the dedicated source. The source, installed lock and CI pin
agree. Offline verification passed for **199 core files, 112 fixtures
and 369 evidence files**.

## M5 monthly TR3 and RSA3 at the R boundary (2026-11-24)

Version `0.1.0.9015` accepts TR3/RSA3 in preprocessing and complete
seasonal adjustment. Six automatically identified models reproduced the JVM
orders, fitted AR/MA and regression coefficients, log/level choice,
likelihoods, **573 dated residuals** and **3,960 values across 60 dated
component rows** bit-for-bit. The fitted-model scope includes the first
reference comparison and controller sweep. A 36-month series verified
that a further identification round raises the core refusal in both APIs.
RSA3 publishes dated residuals; template-specific residual diagnostics
and original-series forecast errors await independent qualification.

Serial and two-worker OpenMP each passed **448/448** RSA3 assertions,
**919/919** RSA1, **130/130** TR1, **498/498** RSA0 and the Java-free
smoke check. All **302,990/302,990** checks passed on the preceding
`0.1.0.9014` serial source suite.
[RSA3 replay](inst/tinytest/test-native-tramoseats-rsa3.R)

## Monthly raw TR3 and RSA3 core integration (2026-11-24)

The package and CI now pin Kalos `0ac66ffaac5c8ff7bb913813485a94da9ed58748`
for first-pass monthly TR3 identification and SEATS extraction from
original observations. The provider builds now include the late controller
and its seasonal F and spectrum sources. Offline verification passed for
**197 core files, 112 fixtures and 368 evidence files**.

## M6 public residual diagnostics for TR1 and RSA1 (2026-11-24)

Version `0.1.0.9014` exposes the preprocessing residual standard error,
the 24-lag Ljung–Box test (22 degrees of freedom) and Doornik–Hansen
normality test (2 degrees of freedom) through both monthly TR1 entry
points. Six original, two shifted level and four monthly TC models
matched all **60 public JVM scalars bit-for-bit**. At 24 residuals,
the Ljung–Box result is undefined and normality remains available.

Serial and two-worker OpenMP each passed all **919/919** RSA1 checks,
**130/130** standalone TR1 checks, **498/498** RSA0 fit checks and
**1,163/1,163** supplied-model checks, plus isolated Java-free smoke
tests. The full serial suite at `0.1.0.9013` passed **302,924/302,924**
checks across 62 files.
[RSA1 diagnostic and component replay](inst/tinytest/test-native-tramoseats-rsa1.R)

## Monthly TR1 residual diagnostics checkpoint (2026-11-24)

The package and CI now pin Kalos
`f65723338f5e7d014d7199d48c4be97385cecb09`. Its public residual
standard error, Ljung–Box and Doornik–Hansen projections match 60
binary64 outputs from twelve independent JVM fits, including the
four monthly TC cases. Offline verification passed for **195 core
files, 112 fixtures and 366 evidence files**.

## M5 monthly transitory changes in R (2026-11-24)

Version `0.1.0.9013` restores stock monthly TC effects in raw RSA1
and caller-supplied TRAMO decompositions. `seats_from_tramo()` accepts
dated TC coefficients with the stock monthly decay of 0.7 on both
scales. The 12 fits match **120 component rows and 7,560 public JVM
values bit-for-bit**, including four cases that estimate TCs.

Serial and two-worker OpenMP each passed **853/853** RSA1 assertions,
all four related numerical suites and isolated Java-free smoke tests.
The complete OpenMP source check at `0.1.0.9011` passed all
**302,491/302,491** assertions with the two known compiler warnings.
[RSA1 product replay](inst/tinytest/test-native-tramoseats-rsa1.R)

## M5 integration checkpoint for monthly transitory changes (2026-11-24)

Kalos `bf5c08f1fafc5c8a6729b94213f40e520e12ef7a` supplies stock
monthly TC restoration on both model scales and a paired four-case
public oracle. Both locks and the CI toolchain agree. Offline verification
passed for **194 core, 112 fixture and 363 evidence files**.

## M5 additive AO and LS effects in R (2026-11-24)

Version `0.1.0.9012` restores additive outlier effects in both
`seats_from_tramo()` and `native_tramoseats_fit(series, "RSA1")`.
Two more raw monthly series select the level transformation and estimate
one AO and an LS/AO pair. The affected series, adjustment, components
and forecasts matched **20 JVM rows and 1,200 values** bit-for-bit.

The complete eight-case RSA1 comparison passed **563/563** assertions
in serial and two-worker OpenMP, including **80 rows and 5,160 values**.
Both builds passed the fitted-model, WK, TR1 preprocessing and RSA0 fit
suites and their isolated Java-free smoke tests. The preceding
`0.1.0.9011` passed **302,487/302,487** full serial checks across 62
files. At `0.1.0.9009` all **301,942/301,942** OpenMP source checks
passed with two existing compiler warnings and zero errors.
[RSA1 product replay](inst/tinytest/test-native-tramoseats-rsa1.R)

## M5 integration checkpoint for additive AO and LS effects (2026-11-24)

Kalos `be98d5405b3f44c30035247779d3b252c1adbb66` restores level-scale
outlier effects in monthly fitted and raw SEATS decompositions. The
source and installed locks agree with the CI toolchain pin. Offline
verification passed for **194 core, 112 fixture and 360 evidence files**.

## M5 monthly RSA1 product fit in R (2026-11-24)

Version `0.1.0.9011` accepts RSA1 and TR1 in `native_tramoseats_fit()`.
The original observations determine log/level selection, airline estimates,
automatic outliers and the dated SEATS product. Six public JVM fits
matched all **60 final component rows and 3,960 values bit-for-bit**.
The custom 3/6 backcast/forecast run matched the separately supplied
mixed AO/LS fitted-model route in the initial, final and merged tables.
The per-template capability inventory identifies the RSA1 dated residual
series and records that preprocessing forecast errors remain unqualified.

Both installed builds passed **419/419** RSA1 product checks,
**130/130** TR1 preprocessing checks and **675/675** RSA0 fit and
diagnostic checks. The fitted bridge and WK suites passed **1,160/1,160**
and **7,676/7,676** respectively. Both isolated Java-free smoke tests
executed the two raw adjustment templates and TR1 preprocessing.
[RSA1 product replay](inst/tinytest/test-native-tramoseats-rsa1.R)

## M5 integration checkpoint for raw monthly RSA1 decomposition (2026-11-23)

Kalos `0e0f2723fef730d5ebbe073309d553e1ab43a765` supplies the public
RSA1 product route and its qualified x86-64 cosine implementation. Both
locks and the CI toolchain pin agree. Offline verification passed for
**194 core, 112 fixture and 358 evidence files**. The earlier
`b6fb5f9` OpenMP source check passed all **301,910/301,910** assertions
with its two existing compiler warnings and zero errors.

## M5 monthly TR1 preprocessing from observations (2026-11-23)

Version `0.1.0.9010` adds `native_tramo_preprocess(series, "TR1")`,
also accepting the template name RSA1. It returns the observed-data
log/level decision, airline coefficients, mean and detected outliers,
likelihood and dated preprocessing residuals. The six JVM fits cover
both transforms, five AO/LS effects and 570 monthly residuals. All
coefficients, likelihood fields and residual values matched bit-for-bit.

Both installed builds passed **130/130** preprocessing assertions,
**675/675** RSA0 fit and diagnostic assertions, **1,160/1,160**
fitted bridge assertions and **7,676/7,676** WK extraction assertions.
Both isolated Java-free smoke runs executed TR1 preprocessing and RSA0
product extraction. The public seasonal-adjustment template remains RSA0.
The preceding `0.1.0.9009` passed **301,938/301,938** serial assertions
across 60 files.
[TR1 reference replay](inst/tinytest/test-native-tramo-preprocess.R)

## Incremental installs after a Kalos pin change (2026-11-23)

The new fitted TRAMO record changed the C++ structure layout. An incremental
install reused an older R wrapper object and crashed on the previously
qualified fitted SEATS entry. Both Makevars templates now rebuild all
wrapper and vendored objects when the generated core revision changes.
Repeating the two incremental installs in their existing build trees
rebuilt the stale wrapper and passed all **1,160/1,160** fitted SEATS
checks plus **7,676/7,676** fixed-model WK checks in each build.

## M5 integration checkpoint for raw monthly TR1 preprocessing (2026-11-23)

Both locks and the CI toolchain now pin Kalos
`33bca9b5bbff3d57a18bf3cb6d75653c288aecad`. The checkpoint supplies
the public six-case TR1 preprocessing reference and original-series
estimator. Offline verification passed for **194 core, 112 fixture and
356 evidence files**. RSA0 remains the delivered raw decomposition template.

## Windows checkout integrity for locked snapshots (2026-11-23)

The Windows CI checkout converted the pinned Kalos text files to CRLF
before `sync-core.py --verify`. A fresh Git checkout with
`core.autocrlf=true` reproduced checksum changes in all 194 locked core
files. Both locked export directories now preserve the exact blob bytes.
A fresh Windows-style checkout passed verification for **194 core, 112
fixture and 355 evidence files** while ordinary package sources still
followed the platform checkout setting. At the previous `b6fb5f9`
checkpoint the complete serial suite passed **301,906/301,906**.
[Windows failure](https://github.com/statistikat/rkalos/actions/runs/45265672229/job/138011798506)

## M6 RSA0 residual normality in R (2026-11-23)

Version `0.1.0.9009` returns the public Doornik–Hansen residual normality
statistic, upper-tail probability and two degrees of freedom in
`fit$residuals$doornikhansen`. On AirPassengers the statistic is
7.770552141622452 and the probability is 0.020542156888630703. All nine
statistic/probability pairs matched the two-process JVM capture bit-for-bit.

Both the installed serial and two-worker OpenMP packages passed
**32/32** normality assertions, **643/643** existing RSA0 fit, residual
and forecast-error assertions, **7,676/7,676** WK and **1,160/1,160**
fitted bridge assertions. Both passed the isolated Java-free smoke test.
[Public normality replay](inst/tinytest/test-native-tramoseats-normality.R)

## M6 integration checkpoint for residual normality (2026-11-23)

Kalos `9f6f74e417ec6f3f88ff4e994152a40438303169` adds the public RSA0
Doornik–Hansen result and paired JVM evidence. Both locks and the CI
toolchain pin agree. Offline verification passed for **194 core, 112
fixture and 355 evidence files**.

## M5 monthly SEATS legacy retries in R (2026-11-23)

Version `0.1.0.9008` extends `seats_extract_wk()` and `seats_from_tramo()`
to the monthly airline seasonal MA replacement and the mean-free
(2,1,0)(1,1,0) seasonal AR-to-MA refit. The fitted mixed AO/LS example
now returns its dated adjusted series, trend, seasonal and irregular
components at both the default 0/24 and custom 3/6 backcast/forecast
horizons. The degenerate airline case returns the public matrix-error
status with empty component tables. The raw-series template remains RSA0.

The installed serial and two-worker OpenMP packages each passed
**7,676/7,676** fixed-model extraction assertions across 36 cases and
eight extra airline threshold requests. The 1,095 frozen component rows
contain **58,020** bit-exact values. All six fitted bridge configurations
passed **1,160/1,160** assertions, including 7,338 WK values and 4,005
merged product values matched bit-for-bit to the JVM. Changing the mean
on the mixed model retains the documented refusal. Both builds passed
the four RSA0 fit, residual and forecast-error suites (**643/643**) and
the isolated Java-free smoke test.
[Extraction replay](inst/tinytest/test-native-seats-extract-wk.R) and
[fitted product replay](inst/tinytest/test-native-seats-from-tramo.R)

The preceding version `0.1.0.9007` passed the full serial suite
(**301,157/301,157** across 59 files). The OpenMP source check of
`0.1.0.9006` passed all **301,120/301,120** assertions with its two
existing macOS compiler warnings and zero errors.

## M5 integration checkpoint for monthly SEATS approximation (2026-11-23)

Kalos `c5bcf978cc7e92d0d5d9992b3c2caa6ad9d0eb65` extends legacy
monthly approximation for two model families. Both locks and the CI
toolchain pin agree. Offline verification passed for **192 core, 112
fixture and 350 evidence files**, including the monthly airline threshold
and seasonal AR-to-MA refit captures.

## M6 RSA0 residual independence diagnostic (2026-11-23)

Version `0.1.0.9007` returns `fit$residuals$lb` with the public
Ljung–Box statistic, upper-tail probability, 22 degrees of freedom and
24 lags. This public residual test uses zero-mean autocorrelations.
The statistic and probability are undefined for 24 or fewer residuals.
The installed serial and two-worker OpenMP packages each passed
**41/41** explicit diagnostic checks
covering the six original series and three short public fits, alongside
**42/42** residual, **62/62** forecast-error and **498/498** fit checks.
All seven finite pairs matched the paired JVM capture bit-for-bit.

At `9e87b84`, the complete serial suite passed **301,116/301,116** over
58 files. A fresh OpenMP run at `e92f076` using independent R processes
also passed **301,078/301,078** across all 57 files in a single run.
[Public diagnostic replay](inst/tinytest/test-native-tramoseats-ljung-box.R)

## M6 integration checkpoint for public residual independence (2026-11-23)

Kalos `aea3d152c63ca365b8dc53887aa16048c38323dd` supplies the public
RSA0 Ljung–Box diagnostic and its nine-case JVM evidence. It also
includes the faster original-series forecast error projection. Offline
verification passed for **192 core, 112 fixture and 346 evidence
files**, with matching source and installed locks and CI toolchain pin.

## M6 RSA0 preprocessing residual diagnostics (2026-11-23)

Version `0.1.0.9006` publishes `fit$residuals$tsres` as a dated monthly
series and `fit$residuals$ser` as its scalar standard error. The six
public JVM series begin 13 months after the input start and end with
the input. All **570 residuals and six standard errors** matched
bit-for-bit. The serial and two-worker OpenMP packages each passed
**42/42** residual assertions, **62/62** forecast-error assertions and
**498/498** raw RSA0 fit assertions. A shifted calendar and changing
forecast/backcast spans also preserved the residual domain.
[Residual replay](inst/tinytest/test-native-tramoseats-residuals.R)

The preceding `e92f076` package passed all **301,074/301,074** serial
assertions across 57 files. OpenMP passed **301,078/301,078** across the
same files, including one isolated rerun after the forked test worker
could not open OpenMP shared memory. Clean source checking reported
the existing two compiler warnings and zero errors. Both isolated
Java-free smoke runs passed.

## M6 integration checkpoint for RSA0 preprocessing residuals (2026-11-23)

Both locks and the CI tools pin Kalos `5c95dc158324e76243d1d7e161d993651edff57b`.
The core exposes dated public RSA0 residuals and their scalar standard
error. Offline verification passed for **192 core, 112 fixture and 345
evidence files**. The checkpoint includes the six-case public JVM
residual capture and intervening live automatic-TRAMO estimation changes.

## M5/M6 RSA0 original-series forecast errors (2026-11-23)

Version `0.1.0.9005` exposes the preprocessing `y_ef` projection as a
dated `series_forecast_stdev`. Its companion status identifies available,
zero-horizon and unavailable results. WK forecast uncertainty continues
to appear in the `final` component table.

The installed serial package passed **62/62** assertions on the paired
six-case JVM capture and custom horizons. All **144** preprocessing
values met the core's established `2e-8 * (1 + abs(reference))`
tolerance. The largest absolute error was **2.42e-13**. The **144**
corresponding WK series forecast standard deviations matched bit-for-bit.
The adjacent raw RSA0 and fixed WK suites passed **498/498** and
**7,369/7,369**.
[Forecast-error replay](inst/tinytest/test-native-tramoseats-forecast-errors.R)

## M5 integration checkpoint for preprocessing forecast uncertainty (2026-11-23)

Both locks and the CI tools use Kalos `2fe3d8e8c8b6bfeae9fa4224517e9e0c12228fe8`.
The core computes RSA0 preprocessing forecast standard deviations separately
from the SEATS component uncertainty. Offline verification passed for
**192 core, 112 fixture and 344 evidence files**, including the six-case
paired JVM forecast-error reference.

## M5 raw-series TRAMO/SEATS RSA0 (2026-11-23)

Version `0.1.0.9004` exposes `native_tramoseats_fit()`. On complete
monthly RSA0/TR0 data, it estimates the mean and airline MA coefficients
and returns the likelihood, dated WK uncertainty, additive final
components and horizon series. Four observed-domain aliases preserve the
input's exact time attributes. Other templates, missing or quarterly
data, insufficient residual degrees of freedom and unsupported native
decompositions remain outside this gate. Product consumers can query
`capabilities()$tramo_seats_fit$templates`.

The serial and OpenMP (two-worker) installed packages each passed
**498/498** RSA0 assertions. All **696 conditional standard deviations**
in the nine observed and forecast uncertainty rows matched the
independent WK reference bit-for-bit for the TR0 AirPassengers fit. All
**30 fitted values and 3,960 product values** matched the six public JVM
cases bit-for-bit. Custom spans, zero-horizon output, date preservation
and refusal paths passed. A clean OpenMP source package check completed
with zero errors and the existing Clang/R-header and floating-point flag
warnings. The package overview example reproduced the full observed
calendar and a 24-month trend forecast. The complete native smoke script
passed in serial and OpenMP libraries isolated from the six forbidden
Java-related R packages. The hosted CI smoke step also executes the full
RSA0 fit before checking which namespaces it loaded. Adjacent fitted-model, WK and first-round
differencing tests passed **718/718**, **7,369/7,369** and **142/142**.
At the preceding published checkpoint `5f57f09`, the complete OpenMP
suite passed **299,658/299,658** across 53 files.
[Fit API](man/native_tramoseats_fit.Rd) · [JVM replay](inst/tinytest/test-native-tramoseats-fit.R)

## M5 integration checkpoint for ordinary monthly RSA0 (2026-11-23)

The checked snapshot advances to Kalos `fa7edeae3c68862e163633c7debce3d102248223`.
It includes the raw-series TR0 estimator, its composed RSA0 product entry
and the frozen six-case x86 JVM oracle. Both locks and the CI toolchain
pin agree. Offline verification passed for **192 core, 112 fixture and
343 evidence files**.

## M5 SEATS from an existing fitted TRAMO model (2026-11-23)

`seats_from_tramo()` accepts a complete monthly series, fitted SARIMA
coefficients, fitted likelihood counts and optional AO/LS coefficients.
It returns dated initial and final WK tables with their uncertainty, plus
the product components after deterministic-effect restoration.

The installed serial package passed **718/718** wrapper assertions.
Four admitted fits reproduced all **8,160 binary64 values** across 116
WK and merged rows in the pinned x86 JVM profile. Both mixed-model
horizons raised the expected approximation refusal. The shared dating
and domain-validation regressions also passed **7,369/7,369**.
Ordinary-data parameter estimation has its own subsequent gate.
[API](man/seats_from_tramo.Rd) · [JVM replay](inst/tinytest/test-native-seats-from-tramo.R)

## M5 integration checkpoint for the fitted SEATS pipeline (2026-11-23)

Kalos `c0594b9018d30b747d8834fab96391ceeb87f901` adds the deterministic
component merge, its paired JVM architecture profiles, and the composed
pipeline for an already fitted TRAMO model. Both installed and source
locks agree. Offline verification passed for **190 core, 112 fixture
and 342 evidence files**. The source CI workflow checks the matching
toolchain revision.

## M5 TRAMO differencing decision (2026-11-23)

`tramo_differencing_decision()` runs the ordinary-input 23-C first-round
module on complete monthly and quarterly series after outlier adjustment.
It accepts the builder's named SARIMA orders and reports the resulting
orders, mean, mean-test statistic, fallback decisions and processing log.
The installed serial package passed **142/142** assertions. All 22
captured direct and battery-entry decisions matched the JVM, including the
binary64 mean-test statistics. Short-series airline and already-airline
branches and invalid input checks also passed. The subsequent
controller battery and complete TRAMO fit have separate gates.

The full serial suite at `6b47f31` passed **299,515/299,515** assertions
across 52 files. A clean OpenMP source build at `5f57f09` completed
`R CMD check --no-tests --no-manual` with the two documented compiler
warnings and zero errors.

## M5 standard TRAMO model preparation (2026-11-23)

`tramo_prepare_model(series, template)` invokes the pinned 23-A builder and
returns the initial log/mean flags, SARIMA orders, registered variable
names, one-based missing positions and interpolated monthly series. It
accepts the seven standard templates without custom overrides. Template
expansion here grants no automatic fit or calendar-selection capability.

The installed serial package passed **139/139** assertions. The frozen JVM
replay covered all seven templates, leading/interior/trailing missing
runs, a non-January start, bit-exact interpolation and the refusal cases.
[Preparation API](man/tramo_prepare_model.Rd) · [Builder regressions](inst/tinytest/test-native-tramo-prepare.R)

## M5 fixed-model SEATS component extraction (2026-11-23)

Version `0.1.0.9003` exposes `seats_extract_wk()`. The R gate compares all 36 frozen 24-C WK cases with the structured 24-D
extraction API at Kalos `ef98f14`. It covers component and information
order, absent versus present-empty rows, monthly dates, raw binary64
values and the core's error and refusal decisions.

The independent serial and OpenMP installations each passed **7,369/7,369**
assertions. The 1,055 structured component rows accounted for 55,812 exact
binary64 values, including conditional standard deviations and forecast and
backcast rows. Each positive-length row carried its monthly `ts` domain.
The run included all six refusal cases and the captured decomposition-error
statuses. The primitive needs a supplied fitted model and positive innovation
variance. End-to-end TRAMO/SEATS fitting has its own gate.
[Installed API](man/seats_extract_wk.Rd) · [WK fixture replay](inst/tinytest/test-native-seats-extract-wk.R)

## M5 integration checkpoint — fitted SEATS bridge (2026-11-23)

The separate generated snapshot advances to Kalos `ef98f14`, with
188 core files, 112 fixtures and 335 evidence files. Offline verification
passed. This revision contains the structured WK decomposition, the
qualified fixed-model bridge, and the paired JVM capture for the component
merge. The CI provider recipe follows the new lock. Binding and product
gates remain separate from this snapshot update.

## Automatic X-13 forecast and backcast tables (2026-11-23)

Version `0.1.0.9002` projects the positive-length 21-H.4 windows from Kalos
`0e064dd` into R. The raw series is in `preadjustment$tables$a1a` (forecast)
or `a1b` (backcast). The corresponding `finals$tables` entries use the roots
`d11`, `d12`, `d16`, `d18` and suffix `a` or `b`. Each series carries
its monthly calendar. Zero-length windows omit the optional keys.
`capabilities()$automatic_fit$x11$horizon_tables` lists the ten table names. Observed outputs preserve the source `ts` time
attributes. The four frozen horizon cases passed all 668 window, value,
refusal and calendar assertions. The 13-case automatic suite passed all
1,213 assertions, including zero-horizon omission and source-calendar
precision. [Horizon fixture replay](inst/tinytest/test-native-x13-horizon.R)

## M5 R checkpoint — TRAMO preprocessing (2026-11-23)

The installed R package exposes `tramo_critical_value()`, `tramo_seasonality()`
and `tramo_transform_decision()` on the existing Kalos `0e064dd` snapshot.
The detector returns both tests and dated differences. The transform decision
returns the model, status, log and level criteria, transform choice and remarks.

The 23-B exit gate passed on macOS ARM64 with R 4.5.2, Clang 23.1.5 and
checksum-pinned LP64 OpenBLAS 0.3.34 in a serial package build. All **196
assertions** passed. The installed-package test replayed 10 critical-value
grid points, five detector series, 10 non-robust decisions and four robust
requests that require the separate AO/LS stage. Captured numerical facts
were compared exactly, including the differenced series, both log/level
criteria and the degenerate Friedman test. Calendar and malformed-input
checks also passed. The test ran without loading Java or JDemetra packages.
An independent OpenMP archive installation on the same host passed all
196 checks with two workers. The source-archive check ran documentation,
compiled-code and example checks (`--no-tests --no-manual`). Escaping the Unicode dash in the X-13
capability string preserved the entire returned list byte-for-byte and
cleared the R source-encoding warning. Two warnings remain from the R
headers' Clang diagnostic and the core-required `-ffp-contract=off` flag.

The source archive built and installed. Offline verification matched all
186 core files, 112 fixtures and 331 evidence files at the existing pin.
The capability list advertises the three primitives. A public TRAMO/SEATS
fit still depends on the 24-D assembly and its end-to-end qualification.
[23-B contract](https://github.com/statistikat/kalos/blob/0e064dd18eb45568fc35e1bed337161598ace36f/docs/contracts/tramo-transform-decision-contract.md) · [R fixture replay](inst/tinytest/test-native-tramo-transform.R)

### Full R suite and exact fixture replay

The serial macOS ARM64 run completed all 51 files. **291,938 of 291,997
assertions passed**. The remaining 59 compare locally regenerated inputs
to the frozen X-11 inputs, in the sigma (51 of 125 cases) and replacement
(8 of 22 cases) tests. All native assertions on the captured inputs passed,
including the exact weights, replacement values, pass tables and X-13
allocation results. The three decomposition-mode suites passed 13,518/13,518.

Before the fixture-reader update, the seven affected files had 5,147 failures.
All seven reproduced with the same counts against the pre-M5 installed
binary. The tests now read the frozen numerical columns with a round-trip
binary64 parser, preserving signed infinities, and run the native routines
on the captured input for each case. The independent regeneration
assertions still compare the local R formulas at full precision.

A separate replay using the pinned glibc sine and cosine kernels matched
all 22 replacement inputs and 123 of 125 sigma inputs. Reading the two
large constants in the sigma generator (`1e300` and `1e298`) at their
correctly rounded binary64 values reproduced the final two cases. The
standard suite continues to report the host-R differences. The fixture
bytes, native calculations and tolerances were preserved.
[Exact fixture reader](inst/tinytest/helpers/fixtures.R) · [Sigma replay](inst/tinytest/test-native-x11-sigma.R) · [Replacement replay](inst/tinytest/test-native-x11-replacement.R)

**Generator portability update.** The installed test harness now supplies the
capture host's trigonometry using the pinned Kalos glibc kernels. It parses
the two large decimal literals in the sigma input generator at their
round-trip binary64 values. All 147 local input regenerations now match the
JVM capture exactly. All **18,717/18,717** assertions in the seven affected
files passed on macOS ARM64, including the 59 previously failing assertions.
The capture scripts and CSV files remained byte-identical to the core lock.
The complete serial suite at `0e064dd` with the horizon and generator
updates passed **292,146/292,146** assertions in 51 files.


## Completed split steps

1. Preserved 102 relevant historical commits and original-to-filtered maps.
2. Added the R-declared nine-routine BLAS/LAPACK provider outside the core snapshot.
3. Adapted every delivered native R/Rcpp entry point through 22b.2 to Kalos.
   Shared window/alignment/regression/Jacobian preparation now runs in C++.
   R retains ts conversion, naming, serialization and result display metadata.
4. Exposed lean public engine interfaces, capabilities and embedded version/provider
   metadata. Required runtime dependencies are only Rcpp and stats, plus base R.
5. Pinned **84 core/provider-notice files** and **112 fixture files** to Kalos
   `4f0625cf44d9cba3c512fdd3d97e1ec4b10132da`. Snapshots are generated and verified
   offline; source installation never fetches GitHub. Source archives include the
   verification script and lock under `inst/core.lock.json`.
6. Migrated 36 native test files and added a boundary/metadata/thread-control test.
   Frozen fixture bytes are untouched; a test-only lexical adapter handles one
   historical Persephone namespace lookup. Small legacy domain/forecast result
   formatters retain serialized compatibility without importing Persephone.

## Qualification

- Source archive `R CMD check`: **266,549 passing assertions / 37 files**,
  zero errors and zero notes. One existing R-header/Clang warning about
  `-Wfixed-enum-extension` remains; it is not a package-code diagnostic.
- All numerical tests and the boundary assertion run without loading Java/rjd3.
- Focused preparation equivalence: **4,281 assertions**.
- Archive audit: every generated core and fixture byte matches its source lock.
  Corruption checks reject changed, missing and unexpected snapshot files.
- R 4.5.2, macOS ARM64, Homebrew Clang 23.1.5, R-configured BLAS/LAPACK,
  two actual OpenMP workers. The same core passes all 37 standalone serial,
  OpenMP and sanitizer programs using OpenBLAS.
- The final source archive reports the configured BLAS library and LAPACK
  version. The serial archive passed **266,545 assertions**, and Persephone
  passed 35 backend assertions with Java-related packages absent.

## Linux (x86_64) baseline qualification — 2026-11-15

First qualification run on a new Linux machine. Nothing here restates macOS
results; only what was executed on this host is recorded.

Environment:

- Ubuntu 24.04.4 LTS, kernel 6.8.0-137-generic, x86_64; g++ 13.3.0
  (Ubuntu 13.3.0-6ubuntu2~24.04.1); R 4.5.3 (x86_64-pc-linux-gnu).
- Installed only into a scratch library (`/tmp/opencode/rlib-b`); no system
  R library was written. Note: this machine's site-wide `Rprofile.site`
  hardcodes `.libPaths()` and ignores `R_LIBS_USER`, so test sessions ran with
  `Rscript --vanilla` and an explicit `R_LIBS` listing the scratch library
  first. Machine configuration, not package code; nothing in the repository
  was changed for it.
- Numerical libraries actually used: R's `-llapack -lblas` resolve through
  Ubuntu alternatives to OpenBLAS 0.3.26 (pthread build); LAPACK is the
  netlib-derived LAPACK 3.12.0 inside `libopenblasp-r0.3.26.so`, reached via
  the `F77_CALL` provider in `src/lapack_r.cpp` (no LAPACKE). Confirmed by
  `core_version()` and `ldd` of the installed shared object.
- `configure`'s C++17 OpenMP compile/link/runtime probe passed (OpenMP 201511
  via libgomp); parallel execution is enabled in this build.

Results:

- `tools/sync-core.py --verify`: passed offline; 84 core files and 112 fixture
  files match commit `4f0625cf44d9cba3c512fdd3d97e1ec4b10132da`.
- `R CMD build` and `R CMD INSTALL` of the archive into the scratch library
  succeeded with zero compiler warnings and zero notes.
- Test suite (`tinytest::run_test_file` over all 37 files, no Java loaded):
  240,563 assertions ran; 239,905 passed; 658 failed; 14 of 37 files aborted
  partway. Every X-11 file, plus inputs, likelihood, linear algebra,
  parallel, prediction-moments, primitives and the package-boundary file ran
  fully green (235,011 assertions); the boundary test confirms rJava/rjd3
  namespaces are not loaded. The 658 failures sit in arma-residuals (288),
  fitted-models (354), identification (2) and optimizer (14). Aborted:
  arma-objective, diffuse, fixed-likelihood, forecast-arma, forecasts,
  information, interpolation, interpolation-arma, joint-information,
  legacy-residuals, regression, regression-information, tramo-information,
  tramo-legacy.
- Not run on this machine: full `R CMD check` of the archive, sanitizer and
  serial-rebuild lanes, the focused preparation-equivalence set, and the
  macOS/Windows lanes.

Blocker (pinned core, not fixed here):

- `src/kalos/src/finite_state_factor.cpp` `checked()` throws
  `underflow_error("finite state factor arithmetic underflows the finite
  numeric range")` when an `long double` intermediate is nonzero but its
  cast to `double` is zero. On x86_64 `long double` is 80-bit extended
  precision; on the macOS ARM64 baseline `long double` is identical to
  `double`, so the same products flush to zero before the guard and the
  evaluation returns `available`. A gdb backtrace for the failing case
  shows `FiniteStateFactor::observation` (finite_state_factor.cpp:135)
  accumulating an exact sum of 6.91e-334 during the diffuse filter of the
  period-12 seasonal AR likelihood (`sarima_innovations` →
  `regression_likelihood_impl`), which casts to 0 and throws.
- Effect on Linux: every likelihood/profile evaluation of a model containing
  seasonal AR/MA factors (periods 4 and 12, airline included) returns
  `numerical_failure`, producing `gradient_unavailable`/`unresolved` fits,
  NULL sensitivities and the aborted/failed test files above. Forcing
  `OMP_NUM_THREADS=1` does not change any of it, so this is a deterministic
  extended-precision semantic difference, not thread nondeterminism and not
  the LAPACK provider. Fixing it requires a Kalos core decision (guard
  semantics or accumulation representation); `src/kalos` is generated and was
  not touched here.
- Secondary divergence: two `test-native-identification.R` assertions fail
  with mean relative difference 1.29e-10 between native `dgesvd` singular
  values of structured seasonal numerator Jacobians and the R `svd()` target
  (tolerance 1e-14). OpenBLAS rounding at borderline conditioning; the
  affected statuses and ranks are correct.
- README example on this machine: `prepare_inputs`, `arma_objective_setup`
  and `fixed_likelihood` behave as documented (`fixed_likelihood(ar=0.3)`
  returns `available`, log-likelihood −680.6269, Q = 4939.992 on the sample
  series), but `arma_optimize` on that level-plus-sinusoid series ends at
  ar = 0.9999972 with status `gradient_unavailable` at the unit-root
  boundary. Not yet cross-checked against the macOS baseline; tracked as an
  open item.

## Linux (x86_64) M1 groundwork — 2026-11-15

Same host and toolchain as the baseline above (Ubuntu 24.04.4, g++ 13.3.0,
R 4.5.3, OpenBLAS 0.3.26 providing LAPACK 3.12.0). The pinned core was **not**
rebuilt or edited here; `src/kalos` is byte-identical to commit
`4f0625cf44d9cba3c512fdd3d97e1ec4b10132da` and the upstream underflow-guard fix
is still pending integration. Only what was executed on this host is recorded.

### Archive `R CMD check`

- `R CMD build .` then `R CMD INSTALL` into a scratch library: clean, **zero
  compiler warnings**. The macOS `-Wfixed-enum-extension` R-header warning does not
  occur with g++ 13.3, so no warning of any kind is left on this platform.
- `R CMD check --no-tests` (default options, `_R_CHECK_FORCE_SUGGESTS_=false`, PDF
  manual built): **Status OK** — zero errors, zero warnings, zero notes; the only
  non-OK line is the informational installed-size report (58.3 Mb, of which 27.5 Mb
  are the shipped test fixtures). Documentation, NAMESPACE, Rd usage, examples,
  portability, `$(BLAS_LIBS)`/`PKG_*FLAGS` and compiled-code checks all pass. **No
  repository-level check issue had to be fixed.**
- Full `R CMD check` including tests: **1 error, 0 warnings, 0 notes**. That single
  error is the known pinned-core blocker: `tests/tinytest.R` stops on the first test
  file with `Error: finite state factor arithmetic underflows the finite numeric
  range`, so the check does not reach a per-file report. Known, not new.
- No check issue introduced by this milestone: the two archived runs differ only by
  the test phase.

### Java independence

`R CMD check --no-tests` was repeated with `JAVA_HOME` unset in a library set where
`rJava`, `rjd3*` and `persephone3` are **not installed** (only `Rcpp`, `tinytest`
and base R): **Status OK**, exit 0. In that same session the boundary, ten X-11 and
five preparation/native test files run green with `package:rJava` absent from the
search path. This machine's default `java` is OpenJDK 21.0.11 (17 and 8 also
installed); neither the check nor the tests consult any JVM.

### Serial (no-OpenMP) lane

- `configure` accepts no arguments and reads no dedicated switch: its single lever
  is the C++17 OpenMP compile/link/runtime probe, which by design honours user
  Makevars overrides. The serial lane was therefore selected the supported way — a
  `R_MAKEVARS_USER` file blanking `SHLIB_OPENMP_CXXFLAGS` — so the probe fails and
  `configure` prints "OpenMP unavailable for this C++17 toolchain; using serial
  execution". A `--disable-openmp`/`KALOS_OPENMP` style flag does not exist; adding
  one is an open usability item.
- Verification that the lane really is serial: `src/Makevars` gets empty OpenMP
  flags, no compile line carries `-fopenmp`, the installed shared object has no
  `libgomp` dependency and no `GOMP_*`/`omp_*` dynamic symbols (the OpenMP build has
  all three), and `core_version()$parallel` reports `compiled = FALSE`,
  `openmp_version = 0`. Installed separately in `rlib-b-serial`.
- Whole suite in the serial lane: **240,565 assertions, 658 failures, 14 of 37 files
  aborting**; the OpenMP lane re-measured here gives **240,563 / 658 / 14**, matching
  the baseline above. Failure and abort sets are identical, and so are the abort
  reasons — the blocker is deterministic extended-precision semantics, not threading.
  The only per-file difference is `test-native-parallel.R` (37 assertions serial,
  35 with OpenMP, from its `compiled = FALSE` branch), zero failures in both.
- X-11 lane plus preparation files (12 X-11 files, plus inputs, primitives,
  likelihood, linear algebra and prediction moments): **234,968 assertions per lane**
  (234,320 X-11 + 648 preparation), identical per-file counts, zero failures.
- Equality stronger than counts: each of those 17 files was re-run under a scratch
  harness that hashes the serialized bytes of every value an assertion compares. All
  17 file digests match between the two lanes across 408,959 recorded values, so
  serial and OpenMP return byte-identical results rather than merely equally-rounded
  ones. (The harness itself mis-scores 4 argument-validation assertions in
  `test-native-prediction-moments.R`, identically in both lanes; the official
  tinytest run reports no failure there.)

### SVD fixture tolerance

`tools/svd-tolerance-evidence.md` (with `tools/svd-tolerance-probe.cpp` and
`tools/svd-tolerance-probe.R`) analyses the two remaining
`test-native-identification.R` failures. Summary: they are line 18
(`expect_equal(result$singular_values, values, tolerance = 1e-14)`) in the
`ar = -0.6, difference = 1e-6` case at period 4 and 12 — proven by relaxing that one
line, which takes the file from 2 failures to 0 — and they concern the small
3x2 near-redundant Jacobians (`kappa_2 = 2.72e6`), not the structured seasonal ones,
whose `kappa_2` is 2.98 and which pass at 1e-14. R's `svd()`, a direct `DGESVD` and a
direct `DGESDD` agree **bitwise** on identical inputs here, and the native values are
bitwise `DGESVD` of the native-normalized Jacobian, so the LAPACK target is excluded.
The gap is 6.73e-17 absolute — 0.21 ULP of `sigma_max`, inside the Weyl bound from the
single 0.8-ULP input difference between native and fixture matrices — amplified by the
conditioning into the reported 1.29e-10 mean *relative* difference, whose theoretical
floor is `kappa * eps = 6.0e-10`. **Recommendation: change the fixture tolerance**
(compare singular values against `64 * .Machine$double.eps * max(values)` in the norm,
or derive the tolerance from `kappa`), **keep `DGESVD` as the target**. The test file
is rkalos-owned (`core.lock.json` generates only `src/kalos/**` and
`inst/tinytest/fixtures/**`), so this needs no Kalos round-trip and touches no frozen
fixture byte. The test was **not** changed in this milestone.

### `arma_optimize` unit-root boundary (README example)

Root-caused on Linux; it is **not** the pinned-core long-double issue. The README fit
ends at coordinate 6.7469984232914353 (`ar = 0.99999723157599829`, value
1.1100490271451684) with status `gradient_unavailable`, 257 evaluations, and a
rejection histogram of `parameter_boundary = 97` — no `numerical_failure` anywhere, no
underflow or overflow guard message in a 121-point scan of the profile, and identical
bytes with 1, 2 and 8 threads and in the serial build. Mechanism:
`arma_optimizer.cpp` searches in `w = tanh(reflection coordinate)`, so the optimizer's
fixed step `difference_step = 1e-4` is applied in a saturating variable where
`dc/dw = 1/(1-w^2) = 1.8e5` at that point; every proposal with `w + h >= 1` is rejected
as `parameter_boundary`, which removes the central stencil, and the surviving
one-sided stencils span `c = 4.6 … 6.4` — far too wide for the declared accuracy
budget. All 11 refinement attempts were reproduced exactly: `error` 1.9e1 … 1.0e4
against budgets of 1.4e0 … 1.4e2, all rejected. The objective itself is smooth and
available there (its own second difference falls off as `h^2` down to 1e-12, and a
stencil taken in the reflection coordinate gives `error = 7.7e-8` against a budget of
`3.9e-5`, i.e. accepted), and the profile minimum is interior: a Newton step from the
stopping point gives `c = 6.7655` (`ar = 0.99999733`), so the optimizer arrived at the
optimum but could not certify its gradient. This is boundary/coordinate handling in the
optimizer wrapper (Kalos side), double-only arithmetic, and expected to reproduce on
macOS; that cross-check still has to be run there.

### Not run on this machine

Sanitizer lanes, the focused preparation-equivalence set as its own target, Windows,
macOS, and any run against the fixed core (still upstream).

## Core integration checkpoint — 2026-11-15

The pending upstream fixes have since landed and the pin advanced to
`cd004b08bcb6cdf1e5b22b2cf34af721ccbaf06c` (via `tools/sync-core.py`, `--verify`
clean; 88 core + 132 fixture files). Two core commits moved in: `f96603f` (the
ARIMA(0,1,1) underflow guard, which removes all fourteen suite aborts and the
658 stale-slope failures) and `cd004b0` (profiled variances published as binary64
quotients of the named binary64 sum). The second was found through this package:
on `f96603f` the full suite ran with a single remaining failure class — 15 of 64
pinned forecasts disagreed in the `mean` attribute at `tolerance = 0` (1 ULP on
the profiled variance) — a scratch-width redivision in the forecast kernels that
the core's standalone suite did not cover. Fixed in Kalos (39/39 per lane, both
lanes), re-synced, re-qualified here.

- SVD comparison: the recommendation above is applied (commit `0696fb1`);
  `test-native-identification.R` is 196/196.
- Whole suite, OpenMP lane: **266,549 assertions, 0 failures, 0 skipped, 0 errors**
  across all 37 files — the first all-green Linux state.
- Whole suite, serial lane (`compiled = FALSE` recompiled): **266,545 / 0**. The
  four-assertion difference from the OpenMP lane is `test-native-tramo-information.R`'s
  six OpenMP-only thread-reproducibility checks versus `test-native-parallel.R`'s
  two serial-specific checks; both lanes report zero failures, errors and skips.
- Qualification ran from an `R CMD build` tarball because `sync-core.py --verify`
  rejects in-tree object files; the serial lane again used the `R_MAKEVARS_USER`
  probe-override and `LD_LIBRARY_PATH` was unset at test time so the pinned
  0.3.34 provider (not the system BLAS) served every LAPACK call.
- pykalos advanced to the same pin; its wheel rebuilt warning-free and smoke-verified
  (`core_version()` reports `cd004b08…`).

## Item 22b.3a/3b at the R boundary — 2026-11-16

- New native entry points `.native_x11_sigma_weights` (input, frequency, start
  period, explicit `legacy`/`accurate` policy, projection counts and exclusion,
  sigma limits, mode) and `.native_x11_replacements` (input, frequency, final
  weights, projection handling, mode), with internal R wrappers
  `native_x11_sigma_weights()`/`native_x11_replacements()` for ts-domain
  handling, 1-based window/exclusion rows and NA-padded neighbor positions.
  The policy argument has **no default**: which sigma rule becomes the product
  default remains an open Kalos maintainer decision, so nothing is promoted.
  Neither entry point is exported yet; the objects were added to the
  `Makevars(.in)/(.win.in)` object lists (the two new core translation units
  were synced but not yet compiled).
- `test-native-x11-sigma.R`: all 125 pinned JDemetra cases reproduce bit-exactly
  through R under the legacy policy (final and first-pass weights, both
  dispersion-window sets, pass counts, exclusion sets, per-row states), with
  recomputation determinism, state/weight consistency and the accurate-policy
  invariants (finite, in-range, masked rows weight one, no legacy cutoff state).
- `test-native-x11-replacement.R`: all 22 pinned cases reproduce bit-exactly
  including the 22b.3a→3b chain (the computed legacy weights equal the fixture
  weights), both projection-masking captures (`mode` 0/1), neighbor provenance
  and counts, period averages, correction factors, the corrected series and the
  apply-corrections identity, fallback counting, validation errors and
  wrapper output shapes.
- Full suite after the addition: **270,653 assertions, 0 failures** across all
  39 files on the OpenMP lane; serial lane **270,649, 0 failures** — the same
  four-assertion conditional difference as before (six OpenMP-only thread
  reproducibility checks versus two serial-specific parallel-info checks);
  both new files contribute identical counts per lane (3,390 sigma + 714
  replacement). The serial lane was installed from a fresh `R CMD build`
  tarball with the `R_MAKEVARS_USER` probe override (0 warnings,
  `compiled = FALSE`, no `omp_*` symbols).

## Core integration checkpoint and B pass at the R boundary — 2026-11-16

- The pin advanced to `16f27859896cc2c2faae69e72e7c24e088f93153` (via
  `tools/sync-core.py`, `--verify` clean offline): one core commit,
  `16f2785` "Qualify the complete fixed additive B pass (22b.4 S1)". The lock
  delta was verified from the diff: 88→90 core files (added
  `src/x11_b_pass.cpp` and `include/kalos/x11_b_pass.hpp`), 132→138 fixture
  files (`native/x11-bpass-cases.R` plus the five `native/x11-bpass/` files);
  only four pre-existing hashes changed (commit, source manifest,
  `CMakeLists.txt`, `revision.hpp`). The snapshot matches the core's own
  source manifest exactly: every `add_library` translation unit,
  `x11_b_pass.cpp` included, is present in `src/kalos`.
- New native entry point `.native_x11_b_pass` (input, frequency, start period,
  explicit `legacy`/`accurate` policy, fixed per-phase initial and final
  seasonal requests, fixed Henderson length, sigma limits) with the internal
  wrapper `native_x11_b_pass()` for ts-domain handling: the 19 pinned member
  tables come back as `ts` on their original calendars (B2..B4g on the
  cropped span via `.native_time_domain`, as in the initial-table wrapper),
  named per-phase filters align by calendar period, and the correction
  history carries 1-based analysis rows with analysis periods plus
  row-basis/available-row provenance. The policy argument has **no default**.
  Neither entry point is exported; `rcpp_x11_b_pass`/`x11_b_pass` were added
  to the `Makevars(.in)/(.win.in)` object lists and the Rcpp registration was
  regenerated (additive-only diff).
  Boundary notes: the 20th public table (B9d) is not a struct member but the
  second round's detection series, so the wrapper exposes 19 member tables
  plus three rounds carrying detection, analysis position, policy, sigma
  weights, replacements, correction factors and replacement counts — the
  final round's sigma weights and correction factors are exactly the pinned
  B17 and B20. The core's validation vocabulary (monthly/quarterly only,
  in-range first period, odd Henderson length ≥ 5, span gates, per-period
  request support) passes through unchanged.
- `test-native-x11-bpass.R`: all 18 pinned fixed additive B-pass cases
  reproduce bit-exactly through R under the legacy policy — 20 tables per
  case on their pinned spans (NaN-aware), correction-round policies, pinned
  analysis positions, replacement counts cross-checked against the B4/B9
  correction tables, round-3 flagged-row counts against B17, the detection
  and B17/B20 correction-history linkages, the additive composition
  identities recomputed in R, determinism, input immutability, the
  accurate-policy runs, the wrapper calendar/named-filter/plain-vector paths
  and the core's validation errors (1,395 assertions).
- Full suite after the addition: **272,048 assertions, 0 failures** across all
  40 files on the OpenMP lane; serial lane **272,044, 0 failures** — the same
  four-assertion conditional lane difference as before. The serial lane was
  installed from a fresh `R CMD build` tarball with the `R_MAKEVARS_USER`
  probe override (`nm` reports no `gomp` symbols).

## Core integration checkpoint and C pass at the R boundary — 2026-11-16

- The pin advanced to `7e8e38dad6d458e5a2827c6feca8dcf3c9326e80` (via
  `tools/sync-core.py`, `--verify` clean offline): one core commit, `7e8e38d`
  "Qualify the complete fixed additive C pass (22b.4 S2)". The lock delta was
  verified from the diff: 90→93 core files (added `src/x11_c_pass.cpp`,
  `include/kalos/x11_c_pass.hpp` and the new shared internal header
  `include/kalos/detail/x11_pinned.hpp`, which the S2 commit factored out of
  the B pass and which the manifest does export into the snapshot), 138→143
  fixture files (the five `native/x11-cpass/` files); six pre-existing hashes
  changed (commit, source manifest, `CMakeLists.txt`, `revision.hpp`,
  `x11_b_pass.hpp`, `x11_b_pass.cpp` — the S2 commit reworked the B-pass
  sources onto the shared header). The snapshot matches the core's own source
  manifest exactly.
- New native entry point `.native_x11_c_pass` (input, frequency, start period,
  explicit `legacy`/`accurate` policy, fixed per-phase initial and final
  seasonal requests, fixed Henderson length, **the required B-pass companion**,
  sigma limits) with the internal wrapper `native_x11_c_pass()`. Design
  decision on the companion: the wrapper takes the caller's finished
  `native_x11_b_pass()` result object and rebuilds the pinned
  `kalos::X11BPass` from it faithfully (NA marks the pinned missing flags,
  1-based analysis rows back to the 0-based core convention) instead of
  silently recomputing the B pass. Recomputing would double the cost and
  could paper over a stale or mismatched B result; passing the caller's object
  through makes the pinned core's own coupling errors fire exactly as in the
  pinned `X11Kernel` composition (completed-status, setting-echo, round- and
  sigma-policy match, B20 coverage, B1-equals-input). The 14 pinned C tables
  come back as `ts` on their original calendars, `C2`/`C4`/`C5a` on the
  cropped span (length n−2·`c2drop`, `c2drop` = frequency/2) and every other
  table on the full span, with the final-stage C13 correction round (analysis
  position, period, C13/C17/C20 linkage) and row-basis/available-row
  provenance. The policy argument has **no default**; neither entry point is
  exported; `rcpp_x11_c_pass`/`x11_c_pass` were added to the
  `Makevars(.in)/(.win.in)` object lists and the Rcpp registration was
  regenerated (additive-only diff).
- `test-native-x11-cpass.R`: all 18 pinned fixed additive C-pass cases
  reproduce bit-exactly through R under the legacy policy — 14 tables per
  case on their pinned spans (NaN-aware), echoed settings and `c2drop`, the
  pinned flagged/corrected row counts, the round policy, position and
  C13/C17/C20 linkages, the additive composition identities recomputed exactly
  in R (C1 = Y − B20, C4 = crop(C1) − C2, C6 = C1 − C5, C9 = C1 − C7,
  C11 = Y − C10, C13 = C11 − C7), determinism, input immutability, the
  accurate-policy composition, the B-coupling rejections (mismatched round
  policy, broken status, altered B1 series, non-covering B companion, malformed
  companion object), the core's validation errors and the wrapper
  calendar/named-filter/plain-vector/missing-companion paths
  (1,124 assertions).
- Full suite after the addition: **273,172 assertions, 0 failures** across all
  41 files on the OpenMP lane; serial lane **273,168, 0 failures** from the
  same `R CMD build` tarball with the `R_MAKEVARS_USER` probe override — the
  same four-assertion conditional lane difference as before.
- Upstream feedback for the Kalos maintainers (recorded here so the milestone
  report carries it; no contract files were touched):
  - (S1) The public B view lists 20 tables but `X11BPass` holds 19 members —
    B9d lives only as the second round's `detection` series (rounds pinned in
    b4d, b9d, b13 order); the contract should state that mapping explicitly.
  - (S1) B17/B20 are bit-identical duplicates of the final round's sigma
    weights and `replacements.correction_factors` — the replacement values
    themselves exist only in the round history — so the contract should say
    plainly that B20 holds correction *factors*, not replacements.
  - (S1) Gate ordering: the quarterly "three annual radii plus one cycle"
    threshold (n ≥ 10) is shadowed by the trend gate (n ≥ trend_length) for
    every quarterly request with `trend_length > 10`, including the pinned
    fixture default of 13, so that message is only reachable with
    `trend_length ≤ 9`; the monthly threshold (n ≥ 30) surfaces naturally.
  - (S2) The cpass fixture's `inputs.csv` is bit-identical to the bpass
    fixture's (verified with `cmp` against the pinned tree) — the same 18
    series and case ids — so the cpass cases replay directly on the bpass
    case list (rkalos reuses `native/x11-bpass-cases.R`); a note in the
    fixture metadata would save integrators from building a second loader.

## Core integration checkpoint, D pass and selection flow at the R boundary — 2026-11-16

- The pin advanced two core commits to `0ad8c59e36351a0a765e98aa0844377555966ec8`
  (via `tools/sync-core.py`, `--verify` clean offline): `886f4e1` "Qualify
  the complete fixed additive D pass (22b.4 S3)" and `0ad8c59` "Qualify the
  pass-level X-11 selection flow (22b.4 S4)"; the fixtures and tests are
  independent, so one reviewed sync covers both. The lock delta was verified
  from the diff: 93→97 core files (added `src/x11_d_pass.cpp`,
  `include/kalos/x11_d_pass.hpp`, `src/x11_selection_pass.cpp` and
  `include/kalos/x11_selection_pass.hpp`), 143→156 fixture files (the six
  `native/x11-dpass/` and seven `native/x11-selection/` files); seven
  pre-existing hashes changed (commit, source manifest, `CMakeLists.txt`,
  `revision.hpp`, `detail/x11_pinned.hpp` — the S3/S4 commits extended the
  shared header — plus `x11_c_pass.hpp`/`x11_c_pass.cpp`, where S3 added the
  `first_period`/`reference_series` echoes). The snapshot matches the core's
  own source manifest exactly.
- New native entry point `.native_x11_d_pass` (input, frequency, start
  period, explicit policy, fixed seasonal requests, Henderson length, the
  required C-pass companion, backcast/forecast horizons; sigma bounds are
  deliberately not taken — the D pass performs no extreme-value analysis)
  with the internal wrapper `native_x11_d_pass()`. The C companion follows
  the established faithful-rebuild/passthrough pattern: the boundary
  rebuilds the pinned `kalos::X11CPass` from the caller's
  `native_x11_c_pass()` result (including the S3 echoes, which is why the
  raw C result now also carries the `reference_series`) and the pinned core
  enforces every coupling check; the C pass is never recomputed. Extended
  calendars are first-class: tables run over the extended series, the
  wrapper derives the dropped original calendar (row span
  `[backcast+1, n-forecast]`, first period `(first_period-1+backcast) %% f
  + 1`) and all four pinned extended fixture cases reproduce the pinned
  original-calendar MSRE diagnostics bit-exactly. The 16 public tables
  (with the pinned `D9_g_bis` alias of `D9bis`) come back as `ts` —
  `D2`/`D4` on the cropped span — and the D9 MSRE table arrives as a
  per-phase `data.frame` with the count-weighted global ratio.
- New native entry point `.native_x11_selection_pass` (requests in the
  seasonal-plan vocabulary `stable`/`s3x1`/`s3x3`/`s3x5`/`s3x9`/`s3x15`/
  `x11default`/`msr`, one per calendar period; `trend_length` 0 for the
  pinned automatic selection; required sigma policy; `jdplus` or
  `x11_five_year_cap` MSR iteration policy; horizons) with the internal
  wrapper `native_x11_selection_pass()`. Decision records come back intact:
  the four trend-stage records (status/length/ratio — the qualified records
  carry a ratio only where the selection engaged, plus the always-computed
  final `iCRatio`), the global-MSR stage with its bit-exact attempt trace,
  and both resolved seasonal plans. Unsupported runs pass the explicit
  `unsupported:<stage>` status straight through with records only — rkalos
  reproduces no pinned degenerate choice. The nested selected-state `b`/`c`/`d`
  keep the fixed-pass boundary shapes, so the records chain back into
  `.native_x11_c_pass`/`.native_x11_d_pass` (tested: bit-exact reproduction
  of the pre-D-final stages, with the all-MSR seas fill only reaching the
  nested `D10bis`, exactly as the pinned bytecode semantics describe).
- `test-native-x11-dpass.R`: all 22 pinned fixed additive D-pass cases
  reproduce bit-exactly through R — 16 tables per case on their pinned
  spans (NaN-aware), the extended inputs verified as hand-extensions of
  their base cases with the pinned original first periods, the per-phase
  MSRE evolutions/counts/ratios and global ratio, the D9 mask rule and its
  count identity against the consumed C round, the nine additive identities
  recomputed in R, determinism, input immutability, the accurate-policy
  composition, the C-coupling rejections (policy, status, settings,
  reference-series equality, coverage, malformed companion), the core's
  validation errors and the wrapper's extended/cropped/original calendars,
  decorated-companion chaining and named-filter paths (1,638 assertions).
- `test-native-x11-selection.R`: all 15 pinned selection cases — the four
  stage records per case against the pinned `selections.csv` (including the
  inline `ratio >= 1 → period+1` contrast at B7 and the deliberate fixed-
  mode "no ratio" convention), the engaged decisions and bit-exact attempt
  traces against `msr.csv` (19-evaluation walk, exhausted and zero-iteration
  windows, `no_history` markers mapped to zero-attempt statuses), the
  `seasFilter` fill, `d9msriter` and B20-correction counts, the full 22-table
  captured selected states of all 13 available cases, the two NaN controls as
  `unsupported:trend-b7` records-only passthroughs, the five-year-cap
  truncation of the 19-evaluation trace with the identical-prefix/divergent-
  `D10bis` invariant, the fixed-chain reproductions, core argument validation
  and the wrapper records/calendar paths (854 assertions).
- Full suite after the additions: **275,664 assertions, 0 failures** across
  all 43 files on the OpenMP lane; serial lane **275,660, 0 failures** from
  the same `R CMD build` tarball with the `R_MAKEVARS_USER` probe override
  (`nm` reports no `gomp` symbols) — the same four-assertion conditional
  lane difference as before.
- Upstream feedback for the Kalos maintainers from these two boundaries (the
  four earlier items above stay verbatim; no contract files were touched):
  - (S3) The S3 commit extended `X11CPass` with the `first_period` and
    `reference_series` echoes that the D-pass coupling consumes — rkalos'
    chained result objects carry those members in the C-pass result shape.
    The contract should state that the pass structs' echo members are part
    of the composition contract: the S2→S3 growth was a schema change for
    any downstream C-result view, not just a new getter.
  - (S3/S4) The selection fixture's `msr.csv` mixes kernel attempt traces
    with `no_history` termination markers and standalone cross-check rows
    for non-engaged cases, and its `final_filter`/`using_default` columns
    restate the final decision on every row; spelling out the row-to-
    `X11MsrAttempt` mapping in the fixture metadata would save the next
    integrator from mis-comparing the attempt traces.
  - (S4) Unsupported selection runs "keep the stage records", but the
    records for stages past the first refused stage are default-constructed
    (`fixed`, length 0, NaN ratio); the contract should say so, since only
    the records up to and including the refused stage carry information.

## Core integration checkpoint, decomposition modes at the R boundary — 2026-11-16

- The pin advanced four core commits to `cce71723be7e4dd4e2e8e6458cfecf9c71b080e0`
  (via `tools/sync-core.py`, `--verify` clean offline): `04b9030` "Version the
  fixed-pass sweep verifiers with repo-relative fixture roots" (touches only
  the core's `tools/` verifiers, which the R-kind snapshot does not carry),
  `db32989` "Qualify the multiplicative decomposition mode (22b.5 M3a)",
  `4646cca` "Qualify the log-additive decomposition mode (22b.5 M3b-i)" and
  `cce7172` "Qualify the pseudo-additive decomposition mode (22b.5 M3b-ii)".
  The lock delta was verified from the diff: 97→98 core files (added
  `include/kalos/x11_decomposition.hpp`) and 156→226 fixture files — 69 files
  in the twelve new `native/x11-{bpass,cpass,dpass,selection}-{mult,log,padd}/`
  directories plus the shared `native/x11-bpass-mult-cases.R` loader; eleven
  pre-existing core hashes changed (`CMakeLists.txt`, `revision.hpp`,
  `detail/x11_pinned.hpp` and the four pass header/source pairs, which gained
  the mode argument and members) plus the commit and manifest hashes. The
  snapshot matches the core's own source manifest exactly.
- Mode at the boundary: a shared `parse_mode`/`mode_name` pair in
  `rcpp_x11_table_support` accepts the four pinned mode names verbatim. The
  B and selection entry points take the mode as a trailing optional argument
  (`NULL`/absent means additive, so every pre-mode positional call stays
  valid); the C and D entry points have no mode argument because the pinned
  kernel runs them on the B pass's context — their TUs now read the mode
  from the caller's companion result echo when rebuilding the pinned struct
  (objects predating the mode echo are treated as additive, which is exactly
  what they are) and echo `mode` plus a `mode_source` note in the result.
  The raw results carry the mode and a mode-aware scope label. The public
  wrappers follow the same shape: `native_x11_b_pass()`/
  `native_x11_selection_pass()` take `mode = "additive"`, while the C and D
  wrappers take `mode = NULL` (inherit from the companion) and reject an
  explicitly stated mode that contradicts the companion's echo at the R
  boundary — the core offers no restatement surface there, so that
  mismatch error is the boundary's own, clearly labelled as such. Nothing
  is exported; the pinned non-positivity rejection ("multiplicative
  decomposition of non positive series") and the selection pass's separate
  completeness rejection propagate verbatim.
- Mode semantics proven at the boundary: the multiplicative and log-additive
  gates fire on zero or negative rows while flagged-missing rows pass the
  fixed-pass gate (the masked fixture cases) but not the selection pass;
  pseudo-additive legally runs zeros and negatives through the whole chain
  (proven by running a zero-and-negative series through B, C and D while the
  multiplicative gate refuses the same series). The pinned `makePositivity`
  repair keeps B7/C7/D7/D12 strictly positive in the mult and log modes
  (asserted structurally on all cases). `B20` corrections are counted
  against the mode's own neutral value — `!= 1` on the multiplicative factor
  scale (mult and pseudo-additive), `!= 0` on the additive/log scale — as
  the fixtures' `b20nz` column requires. The constant-series degenerate
  controls are mode-dependent: `const-fixed` stays `unsupported:trend-b7`
  in all three modes, but `const-msr` becomes available in the
  multiplicative family, where the ratio-based I/C computation avoids the
  NaN the additive and log-additive flows refuse. Additive table identities
  survive recomputation only in the log-additive mode — where the kernel
  runs the plain additive steps on `log(y)` and relations like
  C1 = log(Y) − B20, D6 = D1 − D5 and D13 = D11 − D12 hold bit-exactly —
  while no additive identity holds in the multiplicative or pseudo-additive
  modes (verified empirically before writing the tests); there the bit-exact
  fixture comparison is the qualification, as documented in the test files.
- New test coverage: `test-native-x11-mode-mult.R` (4,461 assertions),
  `test-native-x11-mode-log.R` (4,701) and `test-native-x11-mode-padd.R`
  (4,356), sharing a new `helpers/mode-qualification.R`. Each mode replays
  the complete pinned family through R: all 24 fixed B-pass cases (19 tables
  NaN-bit-exact on their pinned spans, round replacement counts, positivity
  where applicable, determinism, input immutability), all 24 C-pass cases
  (14 tables, round correction count, mode-inheritance echoes, log-domain
  identities in the log mode), all 28 fixed D-pass cases (16 tables, the
  extended inputs verified as hand-extensions with the pinned original first
  periods, per-phase MSRE evolutions/counts/ratios and global ratio — the
  relative form in the multiplicative modes — the D9 mask count identity,
  C-companion coupling rejections) and all 15 selection cases (four stage
  records, the global-MSR attempt traces bit-exact, the captured selected
  states, nested mode echoes), plus the mode validation and non-positivity
  paths and the wrapper inherit/reject/default paths. The four additive
  pass test files pass unchanged: the additive default is bit-identical.
- Full suite after the additions: **289,182 assertions, 0 failures** across
  all 46 files on the OpenMP lane; serial lane **289,178, 0 failures** from
  the same `R CMD build` tarball with the `R_MAKEVARS_USER` probe override
  (`nm` reports no `gomp` symbols) — the same four-assertion conditional
  lane difference as before. No core-caused mismatch appeared at the new pin.
- Upstream feedback for the Kalos maintainers from the mode boundaries (the
  seven earlier items above stay verbatim; no contract files were touched):
  - (M3a) `check_decomposition` discards its pass-label argument and reports
    the same `"multiplicative decomposition of non positive series"` for
    both gated modes, so a downstream boundary cannot tell which pass or
    which mode refused the input; folding the label (and mode) into the
    message would make the verbatim propagation more useful.
  - (M3a–M3b-ii) The `mode` member added to the fixed-pass structs is
    composition-contract surface: the C and D passes have no mode argument
    and inherit the mode solely through `b.mode`/`c.mode`, so any faithful
    downstream rebuild of a companion must carry the echo (rkalos added it
    to its result shapes for exactly this reason). The contracts should say
    so — and `X11SelectionPass` is the only pass result without a mode echo
    member, which would save integrators from remembering it externally.
  - (M3a–M3b-ii) Two fixture semantics are mode-conditional and only
    discoverable by running: `b20nz` counts deviations from the mode's
    neutral B20 value (1 on the multiplicative scale, 0 on the additive/log
    scale) while the column name stays identical, and the constant-series
    `const-msr` control flips from the unsupported degenerate choice to
    available in the multiplicative family. Stating both in the selection
    fixture metadata would save the next integrator from reverse-engineering
    them.

## Core integration checkpoint and the restricted X-13 fit at the R boundary — 2026-11-16

- The pin advanced three core commits to `2b93ab00ab69708b2e5b871f313bed567a04a686`
  (via `tools/sync-core.py`, `--verify` clean offline): `efd4ce0` "Qualify the
  X-13 deterministic-effect allocation staging (22c M3c)", `f84489e` "Deliver
  the restricted fixed-model X-13 fit composition (22c M3-finale)" and
  `2b93ab0` "Qualify the restricted fixed-model X-13 fit on pinned monthly
  fixtures (22c M3-finale)". The lock delta was verified from the diff:
  98→102 core files (added `include/kalos/x13_allocation.hpp`,
  `x13_fit.hpp` and their sources) and 226→243 fixture files (the 8-file
  `native/x13-allocation/` and 9-file `native/x13-fit/` corpora); only two
  pre-existing core hashes changed (`CMakeLists.txt`, `revision.hpp`) plus
  the commit and manifest hashes. The snapshot matches the core's own
  source manifest exactly.
- The allocation seams: `src/rcpp_x13_allocation.cpp` exports
  `.native_x13_preadjust(values, log, variables)` (the fixed-coefficient
  initialStep) and `.native_x13_restore(pre, d10, d11, d12, d13, c17, mode)`
  (the finals) with the pinned allocation vocabulary — one variable list per
  model variable in pinned model order with `class` (outlier, days_related,
  moving_holiday, user), exact `effect` bucket (undefined, series, trend,
  seasonal, seasonally_adjusted, irregular, calendar_effect), one
  `coefficient` per regression `column` vector over the input calendar.
  The restore consumes the caller's finished preadjustment object through
  the established faithful-rebuild/passthrough pattern: the boundary only
  rebuilds the pinned record and the pinned seam enforces the coupling
  (available status, shared calendar). The R wrappers
  `native_x13_preadjust()`/`native_x13_restore()` decorate the stored and
  final tables on the input calendar (ts start/frequency, `row_basis`,
  `available_rows`) and refuse objects that did not come from the boundary.
- The restricted fit: `src/rcpp_x13_fit.cpp` exports `.native_x13_fit()` and
  the `native_x13_fit()` wrapper composes the qualified pieces exactly as the
  core composes them, mirroring the core's fail-closed boundary verbatim and
  never inventing R-side refusals: only fixed seasonal requests ("x11default"
  or "msr" answer `unsupported:x13 fit fixed seasonal requests only`), a
  fixed odd Henderson `trend_length` (the automatic 0 answers
  `unsupported:x13 fit fixed Henderson trend length`), the fixed-path
  horizon domain `backcast = forecast = 0` (anything else, negatives
  included, answers `unsupported:x13 fit requires zero horizons`) and the
  REQUIRED sigma policy: `sigma = NULL` (the default) is passed through
  untouched — no layer applies a default — and the core answers
  `unsupported:x13 fit requires an explicit sigma policy`; it is never an
  R stop. The guards fire in the pinned order (seasonal, trend, horizons,
  sigma), proven by precedence cases. The provenance echo (frequency,
  calendar, log, seasonal requests, trend length, `sigma_set`/`sigma`
  (NULL when unset), sigma cuttings, MSR policy, `mode_requested`/
  `mode_effective` with the updateSpec forcing note, horizons) is complete
  whether or not the run completed; `$preadjustment` (with its stored
  a-tables and `out`) is attached once past the guards, `$selection` (the
  full selection record of the run over `out`, scope
  "restricted_fit_selection_run") once the selection stage ran — the
  selection's own refusals pass through verbatim with the stage records —
  and `$results` (the D10–D13/C17 exactly as finals consumed them, with the
  log-additive results-builder stage note), `$finals`, `$finalSA`,
  `$finalTrend` and `$finalIrregular` only for `status = "available"`.
- Boundary refactor: the selection serialization moved to a shared TU pair
  (`rcpp_x11_selection_records.{hpp,cpp}`, namespace `kalos_x11_rcpp`) whose
  `selection_result()` builder is used verbatim by the pass-level selection
  entry and by the fit, and the X-13 marshalling moved to
  `rcpp_x13_shared.{hpp,cpp}` (namespace `kalos_x13_rcpp`); on the R side the
  selection wrapper's decoration extracted to
  `.x11_decorate_pass()`/`.x11_decorate_selection_record()` so a fit's
  `$selection` reads bit-and-shape-identically to the selection wrapper's
  own result wherever it is chained.
- Capability honesty (the schema decision): `capabilities()` kept every
  schema-1 field with its schema-1 meaning and was extended additively —
  `schema_version` 1→2 marks the extension. **`fit` stays FALSE** and keeps
  meaning UNRESTRICTED fitting (the ARIMA route: estimation, missing rows,
  nonzero horizons, LengthOfPeriod transforms, benchmarking): the restricted
  fit is gated by its own explicit fields, never by `fit`. New top-level
  `native_x13_fit_restricted = TRUE` (the quick gate) and the
  `restricted_fit` block: `schema_version = 1L` (the block's own version),
  `native_x13_fit_restricted = TRUE` (the fit below is qualified end to end
  on the pinned corpus), `fixed_coefficients_only = TRUE` (every variable
  carries its coefficient; nothing is estimated),
  `fixed_seasonal_requests_only = TRUE` (one pinned filter per calendar
  period; x11default/msr are refused), `fixed_trend_length_only = TRUE`
  (odd Henderson >= 5; the automatic 0 is refused),
  `sigma_policy_required = TRUE` ("legacy"/"accurate"; NULL answers the
  core's fail-closed status and no default exists at any layer),
  `zero_horizon_only = TRUE` (backcast = forecast = 0 is the only delivered
  domain; nonzero horizons leave the fixed path through the pinned
  ExactArimaForecasts route), `monthly_only = TRUE` with `frequency = 12L`
  (qualified on the frequency-12 corpus; other frequencies are untested
  here, not rejected — gating is the caller's decision), `modes` (the four
  qualified decompositions after the updateSpec forcing: additive,
  multiplicative, log_additive, pseudo_additive), `refusals` and
  `entry_points`. The `primitives` inventory gained the x13 names plus the
  two qualified-but-unlisted 22b names (`x11_pass_level_selection`,
  `x11_decomposition_modes`); `reason` now names the restricted-fit scope.
  Every existing test only asserts `!fit` and membership, so the extension
  is safe; persephone3 should gate the restricted path on
  `identical(capabilities()$native_x13_fit_restricted, TRUE)` (plus the
  scope fields it cares about) and treat `fit` as unrestricted-fit honesty.
- Tests: `test-native-x13-allocation.R` (239 assertions) drives all 7 pinned
  seam cases through the raw entries — status/type/log/length echoes, the
  a1 series echo, `out` exact (pure algebra in the seam corpus), the twelve
  buckets bit-exact for non-log models and at the pinned 1e-15 relative
  transcendental boundary (`std::exp` vs `Math.exp`, exactly as the core's
  own corpus judges them) under a log model, all nine finals tables exact
  for every effective mode including the forced ones, input immutability,
  the broken-coupling refusals — plus the wrapper decoration and refusals on
  the monthly case. `test-native-x13-fit.R` (636 assertions) drives all 9
  fixed-route cases end to end (complete provenance and forcing echoes, the
  preadjustment at its corpus policy — `out` at 1e-15 under a log model —
  the selection record with named requests and the forced mode, the results
  tables and finals bit-exact outside the log models and at the pinned
  1e-12/1e-9 boundaries inside them, the triplet proven identical to the
  finals tables, the decorated wrapper path on a ts calendar), the two
  arima-route evidence cases refusing verbatim at their pinned horizons with
  the horizons.csv routing facts asserted, the guard-precedence matrix and
  the sigma-required path through both entries without any R stop, the
  byte-identical X13Fast cross-capture of `f_empty`, and the capabilities
  contract above.
- Upstream feedback for the Kalos maintainers from the fit boundary (all
  earlier items above stay verbatim; no contract files were touched):
  - (M3-fin) The two corpora judge the same `tables.csv` `out` column under
    different policies — pure algebra (exact equality) in the allocation
    seam corpus, but downstream of `backTransform`'s exp (1e-15 relative
    under a log model) in the fit corpus, where the captured `out` flows
    into the X-11 run — and the divergence is only visible inside the core's
    own test code. An integrator extending the seam corpus's bit-exact
    assertions to the fit corpus hits a surprise 1-ulp failure on one row.
    Stating the per-corpus comparison policy in the fixture metadata would
    save the next integrator from rediscovering it empirically.
  - (M3-fin) `X13FitResult` echoes the complete X-11 spec but consumes the
    regression side without an echo: after the fact there is no way to audit
    which variable set a fit actually used (count, per-variable column
    counts, the coefficient sum or any other fingerprint). A cheap
    deterministic fingerprint in the provenance would make fit results
    self-describing for downstream pinning without carrying the columns.
- Full suite after the additions: **290,057 assertions, 0 failures** across
  all 48 files on the OpenMP lane; serial lane **290,053, 0 failures** from
  the same `R CMD build` tarball with the `R_MAKEVARS_USER` probe override —
  exactly the 289,182/289,178 M5 baselines plus the 875 new assertions
  (239 allocation + 636 fit), the same four-assertion conditional lane
  difference as before. No core-caused mismatch appeared at the new pin.

## Maintainer sigma-rule default applied at the R interface — 2026-11-17

- Kalos maintainer decision (2026-11-17): the X-11 sigma-rule product
  default is **legacy** — bit-exact against the pinned
  `DefaultExtremeValuesCorrector` / JDemetra behavior. The layering follows
  the agreed split: the core keeps the policy explicit and still refuses a
  missing policy (unchanged core contract, pin stays `2b93ab0`; the wording
  here is rkalos-side until the core documents the decision, which it will
  mirror when it lands), while rkalos applies the default at the R interface.
  This is pure interface-layer work: no sync, no core change.
- Defaulted entry points (the maintainer default applies when the policy is
  omitted or `NULL`; an explicitly stated `"legacy"`/`"accurate"` always
  wins, validated by the core verbatim; the raw `.native_*` entries are
  unchanged and still demand an explicit policy):
  `native_x11_sigma_weights()` (`policy = "legacy"`),
  `native_x11_b_pass()` (`policy = "legacy"`),
  `native_x11_selection_pass()` (`policy = "legacy"`) and
  `native_x13_fit()` (`sigma = NULL` maps to `"legacy"`).
- Inherited entry points (the pinned kernel shares one sigma context across
  the passes, so these inherit rather than default): `native_x11_c_pass()`
  and `native_x11_d_pass()` now take `policy = NULL` meaning INHERIT from
  the B/C companion result's policy echo (the same mechanism as the
  decomposition-mode echo); an explicitly stated policy is a cross-check of
  that echo and a mismatch is refused at the boundary ("inherits its sigma
  weight policy from the …-pass result"), since the core offers no
  restatement surface. Objects predating the companion's policy echo fall
  back to the maintainer default `"legacy"`.
- `capabilities()`: `restricted_fit$sigma_policy_required = TRUE` is
  replaced by `restricted_fit$sigma_policy_default = "legacy"` (the
  restricted_fit block's own `schema_version` bumps 1L→2L to mark the
  substitution); every other field keeps its meaning. The core-side
  explicitness is preserved underneath and documented in the field comment.
- Tests: the wrapper-level "policy must be stated explicitly" refusals in
  the sigma/B/C/D/selection suites become legacy-default proofs — an
  omitted-policy run must be bit-identical to the explicit-legacy run on the
  same pinned fixture (plus the C/D echo-inheritance and mismatch-refusal
  checks) — while the raw-entry refusals stay. In the fit suite the wrapper
  `sigma = NULL` case asserts the default is applied (status `available`,
  bit-identical tables to explicit legacy), the raw entry's verbatim
  `unsupported:x13 fit requires an explicit sigma policy` stays, and the
  core guard precedence with a missing policy (seasonal, trend, horizons,
  sigma — sigma judged last) moved to the raw entry, where `NULL` still
  reaches the core. The capabilities contract test pins the new field set
  and block version.
- Full suite after the change: **290,073 assertions, 0 failures** across all
  48 files on the OpenMP lane; serial lane **290,069, 0 failures** from the
  same `R CMD build` tarball with the `R_MAKEVARS_USER` probe override — the
  M6 totals (290,057/290,053) plus the 16 net new assertions of the reworked
  default/inherence checks (sigma +2, B +1, C +2, D +2, selection 0, fit +9),
  the same four-assertion conditional lane difference as before.

## Evidence fixture root adopted, core pin advanced to e490169 — 2026-11-17

- Infrastructure sync (no R test or interface code changed): the core pin
  advances `2b93ab0 → e490169c95d287c1e66d0b930bf30ddc10a777e8`, spanning the
  contract-doc slice (`27c1318`/`c2a9772`, whose sigma-default resolution the
  section above already mirrors at the interface), the item-21 roadmap
  (`e6ed3a1`), the fixture-root relocation (`d9a5917`: post-split fixtures
  moved to `validation/evidence/`, historical `validation/fixtures/` root
  immutable) and the 21-A estimation-control (`188f4ed`) and Hannan-Rissanen
  (`e490169`) units.
- `tools/sync-core.py` now handles the core's three-section manifest
  (`core` / `fixtures` / `evidence`). Both fixture sections export into the
  SAME `inst/tinytest/fixtures/` tree at their post-prefix paths, so no R
  test path changed; the export refuses if the two sections ever disagree on
  a shared path. Lock **schema 2**: one snapshot entry per section with a
  `section` name and, for fixture sections, the `source_prefix` provenance;
  entries may share a destination directory and `--verify` merges them with
  conflict detection. `--verify` reads schema 1 and 2; schema-2 locks are
  loudly rejected by older tools ("Unsupported source lock version") rather
  than silently misverified — backward-detectability is honest, not silent.
  The `python` kind excludes both fixture sections and still exports exactly
  the 105 core files (106 with the generated `revision.hpp`) — verified into
  a scratch destination with zero fixture/evidence leakage.
- Deltas: core 14 paths (manifest 101→105; `src/kalos` 102→106 files) —
  modified `backends/openblas.cpp`, `CMakeLists.txt`,
  `include/kalos/detail/lapack.hpp`, `include/kalos/linear_algebra.hpp`,
  `revision.hpp` (generated), `x13_allocation.hpp`, `x13_fit.hpp`,
  `src/linear_algebra.cpp`, `src/x11_replacement.cpp`, `src/x11_sigma.cpp`;
  added `estimation_control.{hpp,cpp}` and `hannan_rissanen.{hpp,cpp}`
  (polynomial-root/DGEEV machinery lives inside the existing
  `linear_algebra`/`arma_parameters` surface — no separate roots file).
  Fixtures section 243→112, evidence section +139, and the inst tree grows
  243→251: the 243 relocated files land byte-identically at their old
  inst paths (git sees zero churn) plus the 8 new 21-A/HR corpus files under
  `native/estimation-control{,-hr}/`. Spot-check: one file per fixture/evidence
  directory (84 dirs) — 82 byte-identical against the pre-move
  `2b93ab0:validation/fixtures/…`, the 2 new corpora pristine against
  `e490169`, lock checksums agreeing throughout; 0 mismatches.
- One boundary-shim addition (the only non-R source change): the advanced
  `src/linear_algebra.cpp` calls the newly declared
  `kalos::lapack::dgeev` (companion-matrix roots) and is compiled into the
  shared library, so `src/lapack_r.cpp` gains the matching `F77_CALL(dgeev)`
  definition — link completeness only; no bound code path reaches DGEEV yet.
  The new 21-A sources are NOT added to `Makevars OBJECTS`: the convention is
  compile-only-what's-bound (as with `x11_decomposition.cpp`).
- Orchestrator decision: **no `.native_` bindings for the 21-A intermediate
  estimator units yet** — estimation control and the Hannan-Rissanen
  initializer are internal machinery whose pinned semantics only become
  observable through specification-level outputs; their R-boundary
  qualification lands with the 21-G/H spec-level slices. Recorded for
  evidence-grade awareness: the core's 21-A.2 JVM corpus work found a
  nondeterministic variance-guard case (1e-9-scale guard comparison); any
  future R-side corpus pinning of these units must either dodge or explicitly
  widen that case.
- Full suite after the sync: **290,073 assertions, 0 failures** across all
  48 files on the OpenMP lane; serial lane **290,069, 0 failures** from the
  same `R CMD build` tarball with the `R_MAKEVARS_USER` probe override —
  identical to the M7 totals, as expected for a zero-R-change sync; the 8 new
  fixture files are consumed by the core's own C++ suites, not by R tests.

## Core pin advanced to e19df2c — shared-unit bugfix regression checkpoint — 2026-11-17

- Infrastructure sync (no R or shim code changed): the core pin advances
  `e490169 → e19df2c7b7ded8fa90424a2945d48a2c43e676b9`, spanning the 21-B
  log/level transformation decision (`a942dcb`) and the 21-C outlier scoring
  and forward selection (`e19df2c`). Manifest 105→111 core / 112 fixtures /
  139→150 evidence: `src/kalos` 106→112 files (added
  `transformation_decision.{hpp,cpp}`, `outlier_scoring.hpp`,
  `outlier_scoring.cpp`, `outlier_detail.hpp`, `outlier_scan_engine.cpp`;
  modified `CMakeLists.txt`, append-only `estimation_control.{hpp,cpp}` —
  the new `optimize_sarima_from` warm-start entry — and the generated
  `revision.hpp`); 11 new evidence corpus files under
  `native/transformation-decision/` and `native/outlier-scoring/`
  (inst tree 251→262). The new sources are not in `Makevars OBJECTS`
  (compile-only-what's-bound), so the shared library content is unchanged.
- **Regression verdict (the point of this checkpoint):** the core range is
  purely additive — 6,362 insertions, 0 deletions across every file; the
  already-qualified shared-unit sources (`arma_*`, `linear_algebra`,
  `x11_*`, `x13_*`) are bit-identical between the pins. The three
  transcription bugfixes announced with 21-C — the banded `lcholesky`
  `cnt − 1` update bound, the `Polynomial`/SarimaModel trailing-zero degree
  strip, and the MaLjungBox `calcg(int)` column-*overwrite* (not `-=`)
  recursion — are pinned **inside the new, R-unbound 21-C unit**
  (`outlier_detail.hpp`/`outlier_scan_engine.cpp`, e.g. the
  "Pinned bound: while (s < k + cnt)" loop): the Ansley/MaLjungBox machinery
  there is a separate transcription, not code the fixed-order kernels call.
  The full R suite confirms this observationally: **290,073 assertions,
  0 failures** (48 files, OpenMP lane) and **290,069, 0 failures** (serial
  lane, same `R CMD build` tarball with the `R_MAKEVARS_USER` probe
  override) — bit-identical totals to the M8 checkpoint, with zero
  previously green fixture shifting. **The three bugfixes have zero
  observable effect at the R boundary for existing bindings.**
- Bugfix awareness for the coming 21-G/H R-qualification: those pinned
  semantics become observable only when the 21-C outlier driver is wired;
  the future R corpus pinning must then deliberately include cases that
  discriminate each fixed behavior (fitted MA > 0 for the `calcg` overwrite,
  trailing-zero/stripped-degree polynomials for the degree normalization,
  deep banded updates for the `cnt − 1` bound) — none of the existing R
  fixtures exercise them. Related evidence-grade note from 21-B: the
  discrete `Multiple` alternate-start sign (+/-0.9) shows basin flips and
  flat-ridge stops on four pinned cases (decision bit still pinned),
  flagged upstream for a 21-A replay-fidelity pass.

## Maintainer license decision: DESCRIPTION-only licensing — 2026-11-18

- Kalos maintainer decision (2026-11-17, standing): R packages carry the
  license ONLY as the DESCRIPTION field `License: GPL (>= 3)` — no LICENSE
  files in R package repositories. The top-level `LICENSE`/`LICENSE.md` were
  already removed (`372cdf1`); this step removes what remained: the vendored
  copies `src/kalos/LICENSE` and `src/kalos/LICENSE.md` (the core's own
  license texts, listed in the upstream manifest's core section).
- Mechanism: `tools/sync-core.py` gains a documented maintainer-decision
  exclusion list `R_KIND_EXCLUDED = {"LICENSE", "LICENSE.md"}` applied to the
  R-kind core allowlist — the excluded entries are neither exported nor
  locked, and the export's stale-entry pruning deletes existing copies, so
  the sync tool converges the tree on its own. Lock (schema 2, pin unchanged
  at `e19df2c`) now records 110 core files instead of 112; `--verify` passes.
- Retained, deliberately: `dependencies/OpenBLAS-LICENSE.txt` and
  `dependencies/lapack-netlib-LICENSE.txt` — third-party license notices must
  ship with the vendored code they cover. The **upstream kalos manifest is
  intentionally NOT changed**: kalos is a C++ library that installs its own
  license, and pykalos keeps its own LICENSE within the same decision scope;
  the exclusion is R-package-shaped, so it lives in the R sync tool. The
  python-kind export is unchanged (verified at `e19df2c`: 112
  `vendor/kalos` files, `LICENSE` and `LICENSE.md` included).
- R-side reference audit: no remaining reference to the deleted files outside
  `tools/migration/history-paths.txt`, the immutable historical record
  (left untouched); `.Rbuildignore`'s `^LICENSE\.md$` entry guards the
  top-level and is retained.
- Gate: fresh `R CMD build` tarball (contains zero `src/kalos/LICENSE*`
  entries) installed with `env -u LD_LIBRARY_PATH R CMD INSTALL` into the
  OpenMP lane; full tinytest suite: **290,073 assertions, 0 failures**
  (48 files) — bit-identical to the M9 checkpoint. The serial lane was not
  rerun: this change deletes two inert text files from the vendored tree —
  observationally inert by construction (no compiled or R-visible content),
  so the parallelism axis cannot discriminate it; the single lane confirms
  the expected bit-identical totals.

## Core pin advanced to 2170435 — 21-D differencing decision, R boundary unchanged — 2026-11-18

- Infrastructure sync (no R or shim code changed): the core pin advances
  `e19df2c → 21704352cfc8c6af45e68041e42a9e7eb6d4f44d` (item 21-D: the
  differencing decision unit plus its support cluster — fdlibm-pinned
  transcendentals, the bit-exact `ssq_levmar`/`ssq_deviances` machines, the
  raw-Hannan-Rissanen pipeline extraction, and the `EstimationSolver` seam).
  Manifest core 111→122; `src/kalos` 110→121 files (added
  `detail/fdlibm.hpp`, `detail/hannan_rissanen_pipeline.hpp`,
  `detail/sarima_stabilization.hpp`, `detail/ssq_deviances.hpp`,
  `detail/ssq_levmar.hpp`, `differencing_decision.hpp` and the matching
  `differencing_decision.cpp`, `fdlibm.cpp`, `sarima_stabilization.cpp`,
  `ssq_deviances.cpp`, `ssq_levmar.cpp`); evidence 150→154
  (`native/differencing-decision/`, inst tree 262→266). The
  `R_KIND_EXCLUDED` license behavior still holds (no vendored license texts;
  lock core section = manifest minus the two license entries plus the
  generated `revision.hpp`). New sources stay out of `Makevars OBJECTS`
  (compile-only-what's-bound).
- **Shared-unit diff-shape analysis** (this range MODIFIED files that
  existed at the previous pin — the red-flag question):
  `estimation_control.hpp` is **append-only** (+16/−0, the new
  `EstimationSolver` seam and `spec_d`/`spec_bd` request fields);
  `estimation_control.cpp` (+347/−11) and `hannan_rissanen.cpp` (+178/−116,
  the four transcription-bug fixes) carry real edits but are **21-A units
  that were never bound and never compiled** — the check against both
  `Makevars` OBJECTS lists (45 core sources) is empty, and the include graph
  of every compiled kernel source and Rcpp boundary file touches none of the
  modified or new headers. All 45 compiled core sources and all headers in
  their include graph are bit-identical between the pins: no red flags.
- **Regression verdict: the 21-D slice has zero observable effect at the R
  boundary.** Full suite on the fresh tarball: **290,073 assertions,
  0 failures** (48 files, OpenMP lane) and **290,069, 0 failures** (serial
  lane) — bit-identical to the M9/M10 checkpoints, no fixture moved.
- Forward note for the 21-G/H R qualification: the four pinned-HR
  transcription fixes (Levinson middle-slot temporary, `npi()` max-of-two
  expression, `updatemodel` `btheta` section start, `finalcorrection`
  `BackFilter` layout) move HR numerics toward the pinned JVM — the 21-A
  verify corpus already drifts accordingly (worst 6.2e-4, still 27/27) — and
  the fdlibm `exp`/`log` kernels with the MXCSR FTZ/DAZ guard will make
  bit-level transcendental semantics observable. Future R corpus
  discrimination cases for the HR/differencing surface must be pinned against
  the POST-21-D core and must exercise each fixed behavior; the upstream
  21-D tier note (the `m_i2_144` step0 root-scan pair drifting 3.2e-5/1.31e-4
  inside a banded 2e-4 tier, everything else bit-identical) is the reference
  shape for what a future R tolerance tier may need.

## Core pin advanced to 7200a56 — 21-E initial ARMA identification, R boundary unchanged — 2026-11-18

- Infrastructure sync (no R or shim code changed): the core pin advances
  `2170435 → 7200a56b0997dcdcf1970514ac44030796a8fd46` (item 21-E, the
  initial ARMA identification unit). Manifest core 122→124; `src/kalos`
  121→123 files (added `initial_arma_identification.{hpp,cpp}`; modified
  only `CMakeLists.txt` and the generated `revision.hpp`); evidence
  154→162 (`native/initial-arma-identification/`, 8 corpus files; inst tree
  266→274). The new unit is unbound and outside `Makevars OBJECTS`.
- Structural no-effect: the range is **purely additive** — 7,389
  insertions, 0 deletions across all files; zero modifications to any
  existing `src/` or `include/` file, so the compiled OBJECTS closure and
  its header include-graph are bit-identical to the previous pin
  (red-flag count 0).
- Observational no-effect: full suite on the fresh tarball —
  **290,073 assertions, 0 failures** (48 files, OpenMP lane) and
  **290,069, 0 failures** (serial lane) — bit-identical to the M9–M11
  checkpoints; no fixture moved. The 21-E unit becomes R-relevant only
  through the 21-G/H specification-level qualification.

## Core pin advanced to 0930c9b — 21-F ARMA refinement, R boundary unchanged — 2026-11-18

- Infrastructure sync (no R or shim code changed): the core pin advances
  `7200a56 → 0930c9b153928735b80a0f4150373216544aeca7` (item 21-F, the ARMA
  refinement/comparison unit `refine_arma_model` on the 21-D differencing
  unit and a new 21-E `select_arma_orders` seam). Manifest core 124→126;
  `src/kalos` 123→125 (added `arma_refinement.{hpp,cpp}`); evidence
  162→170 (`native/arma-refinement/`, 8 corpus files; inst tree 274→282).
- Diff-shape classification: besides the additions, the range modifies the
  unbound 21-D/21-E units — `differencing_decision.hpp` (+21/−15, the
  deterministic record-init fix) and `initial_arma_identification.{hpp,cpp}`
  (+53/−4 / +206/−125, the extracted `select_arma_orders` seam;
  `identify_initial_arma` keeps its semantics through the shared select
  engine). All are R-unbound and uncompiled: OBJECTS red-flag count 0, and
  the include graph of every compiled kernel source and boundary file has
  zero hits on the touched headers. Full suite on the fresh tarball
  (built and verified before the lane runs): **290,073 assertions,
  0 failures** (48 files, OpenMP lane) and **290,069, 0 failures** (serial
  lane) — bit-identical to the M9–M12 checkpoints, no fixture moved.
- Awareness note for the 21-G/H qualification: the `differencing_decision.hpp`
  fix makes the step0 `trial_fit` iteration-record flag deterministic
  (default-initialized instead of an indeterminate value), so the upstream
  step0 root-scan tier is now a **six-case list** — `m_i2_144`, `qtr_rw`,
  `m_ar9995`, `m_double_rw`, `m_trend_seasonal`, `m_n231_rw` (largest drifts
  on `m_i2_144`: 3.2e-5 `rmax` / 1.31e-4 `rsmax`; the other five ≤1.5e-6) —
  replacing the previous single-case wording; future R-side tolerance tiers
  should reference the six-case list and the new 21-F 21-E tight/basin/drift
  panel bands (worst dBIC2 7.9e-06, no decision affected). The 21-F unit
  also pins the re-differencing D-polynomial semantics ((1−B^d)(1−B^period)^bd,
  not (1−B)^d at d=2) — relevant when the refinement stage becomes R-bound.

## Core pin advanced to e7af23f — 21-G regression tests and mean rules, R boundary unchanged — 2026-11-19

- Infrastructure sync (no R or shim code changed): the core pin advances
  `0930c9b → e7af23f222315e320ab30338ff8188e74bd444ed` (item 21-G, the
  regression stage tests and mean rules). The range is **purely additive**
  — 2,704 insertions, 0 deletions; the only `src/`+`include/` change is the
  added `regression_variable_tests.{hpp,cpp}` (+154/+308). Manifest core
  126→128; `src/kalos` 125→127; evidence 170→179
  (`native/regression-tests/`, 9 corpus files; inst tree 282→291).
  OBJECTS red-flag count 0 (no existing file touched; the new unit is
  unbound and uncompiled) — the compiled closure is bit-identical to the
  previous pin.
- Full suite on the fresh tarball (built and confirmed before the lane
  runs): **290,073 assertions, 0 failures** (48 files, OpenMP lane) and
  **290,069, 0 failures** (serial lane) — bit-identical to the M9–M13
  checkpoints, no fixture moved.
- Forward note: 21-G's cl-statistics machine and the mean rules become
  R-observable only at the 21-H integration slice (user-variable detection
  is deferred to 21-I); the new evidence corpus (`cl`/`mean`/`rvt`/`vars`/
  `ucov`/`obs` tables) is the pinning reference for that qualification.

## Core pin advanced to b773d27 — 21-H pipeline slice + FP-contraction pin, R boundary unchanged — 2026-11-20

- Infrastructure sync (no R or shim code changed): the core pin advances
  `e7af23f → b773d27a66ef2ee50445a6930c307111efc99ff8` — a **7-commit,
  64-file range** (+24,397/−578): the cxx_std_17 cleanup `3a54b95`, the
  21-H stage-Z contract and slice-1 decision units (`9a65609`, `0c57c57`),
  the final-information producer 21-H.0 (`48342a2`), the automatic-pipeline
  driver 21-H.1 (`a3bbb4e`), the JVM-oracle qualification 21-H.2
  (`b17ee8d`), and the CI root-cause fix (`b773d27`:
  `target_compile_options(kalos PUBLIC -ffp-contract=off)` for GNU/Clang,
  `/fp:precise` for MSVC — the JVM goldens and the vendored fdlibm kernels
  encode unfused IEEE arithmetic). Manifest core 128→140; `src/kalos`
  127→139 files (added the `arma_diagnostics`, `final_estimator`,
  `final_information`, `special_functions` and `x13_automatic` units with
  `mlj_engine.hpp`, `detail/special_functions.hpp`,
  `detail/filter_utility.hpp` and their headers; modified only the
  **uncompiled** `fdlibm.cpp`, `initial_arma_identification.cpp`,
  `outlier_scan_engine.cpp`, `ssq_levmar.cpp` and the vendored
  `CMakeLists.txt`); fixtures 112 unchanged; evidence 179→190
  (`native/final-information/` 5 + `native/x13-automatic/` 6 corpus
  files; inst tree 291→302). Lock schema 2 with both fixture sections;
  `sync-core.py --verify` green.
- **Structural zero-R-boundary-effect proof** (range classification above,
  computed against both `Makevars` OBJECTS lists): the four modified
  vendored sources are 21-C/21-D/21-E units outside the compiled set, and
  the range adds **no header into the compiled include graph**. The 92-file
  compiled closure (45 core sources + transitive headers, scanned also from
  all 37 rkalos boundary TUs) is **byte-identical between the pins** except
  the generated `revision.hpp` stamp — red-flag count 0; the new units are
  unbound and uncompiled (compile-only-what's-bound).
- **Provider (SCRATCH, rebuilt fresh at `b773d27`)**: `env -u
  LD_LIBRARY_PATH` + explicit pinned OpenBLAS 0.3.34 root; OpenMP tree
  ctest **62/62** and reinstalled into `tmp_opencode/kalos-install`
  (installed `revision.hpp` = `b773d27…`); serial tree ctest **62/62**;
  `-ffp-contract=off` in all 126 `compile_commands.json` entries of each
  tree, warning-free, serial tree with zero `-fopenmp`. Core gates at
  `b773d27`: `check-manifest.py` 140/112/190, `verify-migration.py` 406
  (both unchanged by this range's consumer-side effect).
- **FP-flag reach (verified precisely):** PUBLIC propagation was verified
  on the CMake surface — the installed `kalosTargets.cmake` carries
  `INTERFACE_COMPILE_OPTIONS "-ffp-contract=off"` and a fresh
  `find_package(kalos CONFIG)` consumer compiles **its own TU** with the
  flag (consumer `compile_commands.json` verified; binary prints
  `version 0.1.0 revision b773d27`). The R package is **not** a CMake
  consumer: `src/Makevars` never reads the package config and the vendored
  `src/kalos/CMakeLists.txt` is locked but unused by the R build, so
  rkalos' own compile lines are unchanged (grep-verified zero flag hits in
  the lane install logs). The flag therefore guards the CMake consumer
  surface; on this non-FMA baseline it is additionally bit-neutral for
  everything, consistent with the core evidence.
- Full suite on the fresh tarball (built and confirmed before the lane
  runs): **290,073 assertions, 0 failures** (48 files, OpenMP lane,
  18m 08s) and **290,069, 0 failures** (serial lane, 15m 45s — same
  `R_MAKEVARS_USER` probe override, installed into `rlib-b-serial` with
  all four serial markers: probe fails, no `-fopenmp` line, no
  libgomp/`GOMP_*` in the installed `.so`, `core_version()$parallel`
  `compiled = FALSE`/`openmp_version = 0`). **Zero-drift verdict:**
  per-file assertion counts are *identical* to the M9–M14 checkpoints in
  both lanes, and the cross-lane delta stays the documented conditional
  pair (`test-native-parallel.R` 35→37, `test-native-tramo-information.R`
  1580→1574 = net −4). The 21-H slice and the FP-contraction pin have no
  observable effect at the R boundary; no fixture moved, nothing
  re-pinned.

## Core pin b2573fc + M4 automatic-spec surface, OpenBLAS build policy — 2026-11-20

- Infrastructure sync: the core pin advances
  `b773d27 → b2573fc011aee37eced0e97ddfd86d88349887c2` (item 21-H.3, the
  thirteen-case automatic-pipeline oracle closure) — one commit, 16 files,
  +930/−67. The only `src/`+`include/` changes are `outlier_scoring.{hpp,cpp}`
  (selection-window facts for the probes) and `x13_automatic.cpp`
  (decision-chain recording); all three sit **outside the restricted compiled
  closure**. Manifest counts unchanged (core 139, fixtures 112, evidence 190);
  evidence content updated at the new probe (`oracle-facts.csv` +43 rows of
  21-H.3 fact families, `metadata.R` probe hash `088eb815…`, `SHA256SUMS.csv`
  regenerated). `sync-core.py --verify` green; the vendored provider
  (`src/kalos/backends/openblas.cpp`) and every synced file diff byte-identical
  against `git show b2573fc:…`.
- **Structural zero-drift proof** for the restricted paths (pin diff × the 92-file
  compiled closure × OBJECTS, scanned from all boundary TUs): restricted-closure
  bytes drifting between the pins **NONE** (stamp apart); the three modified core
  files are uncompiled/unbound; the automatic closure adds 39 core sources over
  the restricted one; `src/RcppExports.cpp` +60 and `R/RcppExports.R` +4 are
  additive-only. Red-flag count **0**.
- **Qualification of the vendored pin itself** (fresh CMake tree extracted at
  exactly `b2573fc`, pinned OpenBLAS 0.3.34 provider, Release): **ctest 62/62
  OpenMP and 62/62 serial**, including `test_x13_automatic_oracle` (13 cases,
  4,637 checks, 0 failures). The contract tier bands *are* the b2573fc
  behaviour: the oracle's own worst deviations — diff t-mean 7.42e-4, ARMA
  panel BIC2 8.09e-3, landing coefficients 9.75e-3 — are reproduced to the
  digit by the R lane (below).
- **M4 surface (first R-observable automatic pipeline).** Boundary:
  `native_x13_fit_automatic()` — namespace-internal, **not** NAMESPACE-exported,
  49 arguments; echoes provenance, the module state exactly as received
  (1-based span kept 1-based) and the full decision chain (1-based positions,
  types `AO`/`LS`/`TC`/`SO`); the pinned fail-closed matrix arrives as
  structured `status`, never an error path, never a silent restricted fallback.
  `R/native_x13_automatic_fit.R` decodes the automatic spec with the pinned
  `mirror_request` scalars (max_d 2, max_bd 1, precision 1e-7/1e-5, models 5,
  grid 2/2/1/1, eps 1e-7, max_outliers 30, va 0, check_mu on, accept_airline
  off, ur_limit .95, mean tests 1.96/1.0, tsig 1.0); `ub1`/`ub2` pass as
  reciprocals (fixtures store JVM-space values). `capabilities()` advances to
  schema_version **2** with the `automatic_fit` block (transform/automodel/
  outliers/ami/regression/deferred/x11 surfaces and their required module
  state); `fit` stays **FALSE** until the boundary is user-facing.
- **Tests**: `inst/tinytest/test-native-x13-automatic.R`, 1,167 assertions,
  thirteen corpus cases; all DECISIONS asserted **bit-identical** against the
  pinned JVM goldens (`sigma="accurate"`); tier bands per contract with **no
  relaxation** — measured R worsts equal the C++ oracle measurements: display
  AICC 3.42e-6 (band 1e-5), diff t-mean 7.42e-4 (1e-3), ARMA panel 8.09e-3
  (1e-2), landing coefficients 3.05e-6 (2e-3), ARMA landing 9.75e-3 (1e-2),
  tables 5.04e-7 rel (1e-6), landing stats 4.48e-7 (2e-3).
- **Build policy — the R-vs-C++ tier divergence, root-caused and closed.**
  The `.so` had been binding `cblas_dgemm` (the provider's only dynamic BLAS
  reference) to the *system* OpenBLAS, because R's `$(BLAS_LIBS)` preceded the
  pinned archive at link time: a shared library satisfies an undefined
  reference before archives are scanned, so mixed kernels ran and the
  automatic pipeline became environment-dependent (thread/BLAS knobs moved
  the diff t-mean to 5.6e-3 and ARMA landing to 6.2e-2). R's `$(BLAS_LIBS)`/
  `$(LAPACK_LIBS)` are now **removed** from `PKG_LIBS`; `configure` fails
  closed unless the OpenBLAS C API (`cblas.h`/`lapacke.h` + linkable
  CBLAS/LAPACKE) is reachable — `KALOS_OPENBLAS_ROOT` or
  `KALOS_OPENBLAS_CPPFLAGS`/`KALOS_OPENBLAS_LINK` direct it — and ELF builds
  add `-Wl,--exclude-libs,ALL -Wl,-Bsymbolic` so the object is self-bound
  (no in-process interposition by R's own OpenBLAS). The automatic fit pins
  BLAS kernels single-threaded for the duration of the qualified call (the
  oracle surface is measured per-kernel serial; outer OpenMP untouched).
  After the fix the R lane is bit-identical to the C++ oracle lane **in the
  default environment**, with no thread or BLAS variables set.
- Full suites (fresh tarball lanes, default environment): **291,240
  assertions, 0 failures** (49 files, OpenMP lane — existing 48 files
  per-file-identical to the 290,073 checkpoint, new file +1,167; the
  boundary file's provider string now asserts the vendored C-API provider
  `OpenBLAS`, count 8 unchanged) and **291,236, 0 failures** (serial lane,
  290,069 + 1,167; all four serial markers: probe falls back, zero
  `-fopenmp` lines, no `GOMP_*` in the installed `.so`,
  `core_version()$parallel$compiled = FALSE`).

## Core pin f500162 — 21-H.4 horizon route pinned, fail-closed matrix mirrored — 2026-11-21

- Infrastructure sync: the core pin advances
  `b2573fc → f500162a3be1539803566f9bc8f071e0deb4f66b` (item 21-H.4, the
  nonzero-horizon initialStep route, via the workflow-only 7564ce7) — 26
  files, +3,102/−90 upstream. The six `src/`+`include/` changes
  (`x13_automatic.{cpp,hpp}` horizon wiring incl. the estimation-block AR
  sign correction, `x13_allocation.{cpp,hpp}` allocation-seam extensions,
  NEW `kalos::arima_forecasts`) all sit **inside the compiled closure**
  (`x13_automatic.cpp` has been compiled since M4); `kalos/src/
  arima_forecasts.o` joins OBJECTS — `x13_allocation.cpp` now calls the
  unit (link edge verified). Manifest advances to core **141**, fixtures
  **112**, evidence **195** (five `native/x13-horizon` capture files);
  lock schema 2 unchanged, `sync-core.py --verify` green pre/post, every
  vendored file byte-identical against `git show f500162:…`.
- **Qualification of the vendored pin itself** (fresh `git archive` trees
  at exactly `f500162`, pinned OpenBLAS 0.3.34 provider, Release, built
  with `LD_LIBRARY_PATH` unset, tested with it pointing at the provider's
  `lib/`): **ctest 64/64 OpenMP and 64/64 serial** — incl.
  `test_x13_automatic_horizon` (4 captures, 6,459 checks) and
  `test_arima_forecasts`; the automatic sweep (71,401 checks) and the
  13-case oracle (4,637 checks) pass unchanged, so zero-horizon semantics
  are provably invariant at the new pin.
- **Structural compiled-closure diff b2573fc → f500162** (roots = 64 core
  OBJECTS sources + 37 boundary TUs, transitive include closure, 133
  files): the drifting bytes are EXACTLY the six vendored files —
  `x13_automatic.cpp` (horizon gating + the a1 sign fix; header change
  comment-only), `x13_allocation.cpp` (horizon plumbing, new
  `arima_forecasts` call), `x13_automatic.hpp` (comment-only),
  `x13_allocation.hpp` (extended input structs), +`arima_forecasts.{cpp,
  hpp}` — with **zero unexplained drift**; OBJECTS changed only by the
  additive `arima_forecasts.o`; every boundary TU and both RcppExports
  files are byte-identical to 0a56d1f. Results are nevertheless unchanged
  at the existing surface because the corpus drives `backcast = forecast =
  0` (zero-horizon code paths, invariant per the core's unchanged sweep +
  oracle) and the a1 sign correction is bit-exact per the core's horizon
  and zero-horizon oracles — confirmed at the R boundary below.
- **Pinned fail-closed matrix drift, mirrored with maintainer sign-off.**
  21-H.4 *delivers* the nonzero-horizon route on the automatic driver and
  removed its `unsupported:x13 fit requires zero horizons` gate (the fixed
  fit keeps refusing); the core mirrored its own `t01_status_matrix` in
  the same commit. Three assertion sites in
  `inst/tinytest/test-native-x13-automatic.R` were updated to the new
  pinned matrix, count-preserving (file stays at 1,167 assertions):
  horizon + missing policy now answers the sigma-gate refusal; the
  guard-order case with a horizon present answers the 21-I defer refusal;
  the wrapper horizon case now pins the *shape* of the delivered route
  (`available`, fitted, model/chain non-NULL, horizon echoed in the
  request), not its numbers. `capabilities()` keeps
  `$automatic_fit$zero_horizon_only = TRUE` — R-side horizon qualification
  against the four JVM captures is the next slice; the flag comment now
  states that the refusal itself is gone.
- **Harness + full suites (fresh tarball lanes built from the pinned
  commit, default environment, provider via `KALOS_OPENBLAS_ROOT`,
  object self-bound — zero dynamic BLAS/LAPACK references; R 4.5.3
  toolchain).** 13-case automatic harness: **0 failures**, decisions
  bit-identical to the JVM goldens, tier worsts equal to the b2573fc
  checkpoint to the digit — display AICC 3.42e-6 (1e-5), diff t-mean
  7.42e-4 (1e-3), ARMA panel BIC2 8.09e-3 (1e-2), landing coefficients
  3.05e-6 (2e-3), ARMA landing 9.75e-3 (1e-2), tables 5.04e-7 rel (1e-6),
  landing stats 4.48e-7 (2e-3). Full suites: **291,240 assertions, 0
  failures** (OpenMP lane) and **291,236, 0 failures** (serial lane), both
  **per-file identical to the 0a56d1f checkpoint** (49/49 files). All
  four serial markers hold (probe falls back, zero `-fopenmp` lines, no
  `GOMP_*`/`omp_*` symbols in the installed `.so`,
  `core_version()$parallel$compiled = FALSE`). CI workflows moved to
  `actions/checkout@v5`/`actions/setup-python@v6` (Node 24 runtimes,
  mirroring kalos 7564ce7).

## Item 21-H.4 horizon qualification at the R boundary — 2026-11-21

Core pin unchanged (f500162, no re-vendor): the pinned core already delivers the
initialStep horizon route; this slice qualifies the horizon envelope AT the R
boundary and makes the boundary's decorations horizon-aware.

- `capabilities()$automatic_fit`: schema_version 2L → 3L; the x11
  `zero_horizon_only` marker is REPLACED by `horizon_envelope` (month counts,
  negative = |h| annual frequencies via the core's normalization, downstream
  composition only); the `qualification` string names the new horizon capture.
  `restricted_fit` is byte-identical: the restricted route keeps its {0, 0}
  domain (the core still refuses nonzero fixed-path horizons).
- `native_x13_fit_automatic()`: docstring now states the qualified envelope
  (no more "{0, 0} only" claim). Decoration fixes the horizons exposed: the
  extended composition tables (`out`, `d10..c17`, length `n+bc+fc`) get the
  calendar start shifted by the backcast (mirroring the core's negative-count
  normalization locally), and `row_basis`/`available_rows` follow each table's
  own length instead of pinning every preadjustment table to `n` (which made
  `available_rows$out` overrun into NA indices once a backcast existed). At
  zero horizons both changes are byte-identical to the previous decoration
  (same `ts()` calls, same bases); the restricted wrapper and
  `.x13_decorate_tables` are untouched.
- `test-native-x13-automatic.R` (+20, 1167 → 1187): capability envelope re-pin
  (schema 3L, retired marker, envelope wording), the horizon block upgraded
  from shape-only to envelope acceptance (lengths, shifted tsp starts, row
  basis, available rows, negative-horizon `forecast=-1` → 108 months with the
  raw echo kept), and the fail-closed matrix untouched (horizons are no gate).
- NEW `test-native-x13-horizon.R` (458 assertions): the 4-case JVM horizon
  capture replay (quick_h6, ao_b3, multi_b3h6, so_h6). Evidence gate first
  (SHA256SUMS.csv over the vendored bytes via `tools::sha256sum`, file set
  pinned); then per case: decisions BIT (status, log/mean/orders,
  mode_effective, horizon echo, variable identity, narma, neff), `a1` and
  `d18` BIT claims, and every other boundary-exposed table against the
  kCases band matrix transcribed from `tests/jvm_x13_automatic_horizon.h`
  at f500162 — copied, never widened, gated per case (no cross-case pooling).
  Landing coefficients/loglik/aicc reuse the existing zero-horizon bands
  (2e-3 / 1e-2 / 2e-3, horizon-invariant upstream). The `a1a`/`a1b` and
  `d11a..d18b` slices and `d8` do not cross this surface and stay inside the
  core claim. Measured worsts (all inside band, largest occupancies): out
  1.84e-7/2e-6, 2.69e-7/5e-6, 5.04e-7/1e-5, 1.05e-6/2e-5; d13 2.01e-4 vs
  5e-3 (quick_h6); c17 6.97e-4 vs 1e-2 (so_h6); e2 5.61e-3 vs 1e-1 (ao_b3)
  and 1.65e-2 vs 2e-1 (multi_b3h6); landing worst 3.05e-6 (2e-3).

Shared gates on the edited tree: full suites **291,718 assertions / 0
failures** (OpenMP lane, 50 files) and **291,714 / 0** (serial lane);
**every pre-existing file per-file identical to the f500162 checkpoint**
except the intended `test-native-x13-automatic.R` (+20) and the new horizon
file. The 13-case harness reproduces its recorded worsts to the digit on the
edited lane (3.42e-6, 7.42e-4, 8.09e-3, 3.05e-6, 9.75e-3, 5.04e-7, 4.48e-7;
measured at `OMP_NUM_THREADS=4` as all lane baselines were — an unset thread
count changes only the sub-band tail drift, never a verdict). All four serial
markers hold (probe falls back, zero `-fopenmp` lines, no `GOMP_*`/`omp_*`
symbols, `core_version()$parallel$compiled = FALSE`). Isolated no-Java
verification passes with `--native` on the new build; `R CMD build` clean;
`DESCRIPTION` License remains the sole `GPL (>= 3)`.

## Compatibility

Accurate final-fit covariance and accurate residual/interpolation reporting stay
at their existing defaults; explicit legacy modes remain Java-free. JDemetra MSR
iteration behavior stays the default, with the published cap also available.
Historical `persephone_native_*`/forecast result classes are retained. This
experimental package does not implement full seasonal-adjustment fitting.

## Next work

1. Track integration and platform checks after the initial GitHub publication.
2. Add R conversions/tests only as Kalos delivers each next numerical slice.
3. Keep pinned source/fixture refreshes reviewable, regenerate Rcpp registrations
   when signatures change, and run archive checks and serial/no-Java lanes.
4. Expand platform qualification through the configured Linux/macOS/Windows CI.
   Remote workflow configuration is not a claim that those runs have passed.

The shipped X-13 routes use the documented legacy sigma-weight default at
the R boundary. The next full-method capability gate is TRAMO/SEATS 24-D.
