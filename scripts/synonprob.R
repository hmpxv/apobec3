# observed APOBEC3-type mutations across tree:
obs_syn <- 226
obs_non_syn <- 245
obs_nonsense <- 30
obs_non_coding <- 132
obs_coding <- obs_syn + obs_non_syn + obs_nonsense
obs_total <- obs_coding + obs_non_coding

# contexts in reference genome by the effects of mutation:
exp_syn <- 4990
exp_non_syn <- 14618
exp_nonsense <- 692
exp_non_coding <- 3418
exp_coding <- exp_syn + exp_non_syn + exp_nonsense
exp_total <- exp_coding + exp_non_coding

# probability of hitting a synonymous context
prob_success <- exp_syn / exp_total
successes <- obs_syn
trials <- obs_total

P <- pbinom(successes, size=trials, prob=prob_success, lower.tail=FALSE)
print(paste("P(syn >= ", obs_syn, "): ", P))

# probability of hitting a non-coding context
successes <- obs_non_coding
trials <- obs_total
prob_success <- exp_non_coding / exp_total

P <- pbinom(successes, size=trials, prob=prob_success, lower.tail=FALSE)
print(paste("P(non_coding >= ", obs_non_coding, "): ", P))
