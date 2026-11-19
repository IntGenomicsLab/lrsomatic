// Every process conf/modules.config gives the retry-then-ignore policy, named as in the pipeline
include { PROBE as PADFOOT_SEVERUS_WAKHAN     } from './probe'
include { PROBE as PADFOOT_SAVANA             } from './probe'
include { PROBE as RECONPLOT_SEVERUS_ASCAT    } from './probe'
include { PROBE as RECONPLOT_SEVERUS_WAKHAN   } from './probe'
include { PROBE as RECONPLOT_SAVANA           } from './probe'
include { PROBE as SHATTERSEEK_SEVERUS_ASCAT  } from './probe'
include { PROBE as SHATTERSEEK_SEVERUS_WAKHAN } from './probe'
include { PROBE as SHATTERSEEK_SAVANA         } from './probe'
include { PROBE as CORAL_RECONSTRUCT          } from './probe'
include { PROBE as CORAL_CYCLE                } from './probe'
include { PROBE as CORAL_PLOT                 } from './probe'
include { PROBE as AMPLICONCLASSIFIER         } from './probe'

workflow ERROR_POLICY_PROBE {
    take:
    ch_cases // [meta]: meta.exit = exit status of the first attempt, meta.always = fail every attempt

    main:
    PADFOOT_SEVERUS_WAKHAN(ch_cases)
    PADFOOT_SAVANA(ch_cases)
    RECONPLOT_SEVERUS_ASCAT(ch_cases)
    RECONPLOT_SEVERUS_WAKHAN(ch_cases)
    RECONPLOT_SAVANA(ch_cases)
    SHATTERSEEK_SEVERUS_ASCAT(ch_cases)
    SHATTERSEEK_SEVERUS_WAKHAN(ch_cases)
    SHATTERSEEK_SAVANA(ch_cases)
    CORAL_RECONSTRUCT(ch_cases)
    CORAL_CYCLE(ch_cases)
    CORAL_PLOT(ch_cases)
    AMPLICONCLASSIFIER(ch_cases)

    emit:
    done = PADFOOT_SEVERUS_WAKHAN.out.done.mix(
        PADFOOT_SAVANA.out.done,
        RECONPLOT_SEVERUS_ASCAT.out.done,
        RECONPLOT_SEVERUS_WAKHAN.out.done,
        RECONPLOT_SAVANA.out.done,
        SHATTERSEEK_SEVERUS_ASCAT.out.done,
        SHATTERSEEK_SEVERUS_WAKHAN.out.done,
        SHATTERSEEK_SAVANA.out.done,
        CORAL_RECONSTRUCT.out.done,
        CORAL_CYCLE.out.done,
        CORAL_PLOT.out.done,
        AMPLICONCLASSIFIER.out.done
    )
}
