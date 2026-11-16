source("helpers/fixtures.R", local = TRUE)
source("helpers/mode-qualification.R", local = TRUE)
# Item 22b.5 M3b-i: the log-additive decomposition mode at the R boundary.
# The pinned kernel transforms the input with Math.log before the plain
# additive steps, so the delivered tables are log-domain and the additive
# relations survive bit-exactly on them (C1 = log(Y) - B20, D6 = D1 - D5,
# the D11/D13 final decomposition and the cropped-span differences); the
# positivity gate and the makePositivity repair apply as in the
# multiplicative mode.
local({
  x11_mode_qualification('log','log_additive',
    b20_correct=0,pinned_label='LogAdditive',
    unsupported=c('const-fixed','const-msr'),repair=TRUE,identities='log')
})
