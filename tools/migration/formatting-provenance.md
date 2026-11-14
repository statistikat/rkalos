# Serialization compatibility helpers

`R/result_domain.R` and `R/forecast_contract.R` preserve the plain-list conversion
helpers from Persephone commit `8450d6e`, files `R/result_schema.R` and
`R/forecasts.R`. Functions are scoped with native prefixes in rkalos. They retain
historical result classes/schema and require only base R and stats. They do not
import Persephone, its R6 classes, dashboards, forecasts package or Java adapter.
The original files remain available in the source repository at that commit.
