source("helpers/fixtures.R", local = TRUE)
source("helpers/mode-qualification.R", local = TRUE)
# Item 22b.5 M3a: the multiplicative decomposition mode at the R boundary.
# The full pinned mult fixture family (24 B/C cases, 28 D cases including
# the extended calendars, 15 selection cases) replays bit-exactly through
# the mode-aware pass entry points; the mode is selected on the B pass and
# inherited by the C and D passes, and the pinned makePositivity repair of
# the trend tables keeps B7/C7/D7/D12 positive.
local({
  x11_mode_qualification('mult','multiplicative',
    b20_correct=1,pinned_label='Multiplicative',
    unsupported='const-fixed',repair=TRUE,identities='none')
})
