source("helpers/fixtures.R", local = TRUE)
source("helpers/mode-qualification.R", local = TRUE)
# Item 22b.5 M3b-ii: the pseudo-additive decomposition mode at the R
# boundary. The pinned X11{B,C,D}StepPseudoAdd pipeline runs multiplicative
# factor arithmetic around additive-scale adjustments, so no additive table
# identity survives recomputation (the bit-exact fixture comparison is the
# qualification) and zeros and negative rows are legal inputs: the pinned
# positivity gate is multiplicative/log-only, proven here through R by
# running a zero-and-negative series through the whole chain while the same
# series is refused by the multiplicative gate.
local({
  x11_mode_qualification('padd','pseudo_additive',
    b20_correct=1,pinned_label='PseudoAdditive',
    unsupported='const-fixed',repair=FALSE,identities='none')
})
