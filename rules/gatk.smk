def gatk_extra(picard = False, jdk = False, hc = False):
    extra = ''
    if picard:
        extra += ' --TMP_DIR %s' % config['tmpdir'] #' TMP_DIR=%s' % config['tmpdir']
    else:
        extra += ' --tmp-dir %s' % config['tmpdir']
    if jdk:
        if picard:
            extra += ' --USE_JDK_DEFLATER --USE_JDK_INFLATER' #' USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
        else:
            extra += ' --use-jdk-deflater --use-jdk-inflater'
    if hc:
        extra += ' -pairHMM LOGLESS_CACHING'
    return extra

