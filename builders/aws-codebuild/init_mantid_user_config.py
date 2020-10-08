def _createDefaultMantidUserConfig(facility='SNS', instrument='ARCS', CORES=None):
    # create default Mantid user configuration for DEMO purpose.
    import os
    dotmantid = os.path.expanduser('~/.mantid')
    if not os.path.exists(dotmantid):
        cmd = "git clone https://github.com/yxqd/dotmantid %r" % dotmantid
        if os.system(cmd):
            raise RuntimeError("%s failed" % cmd)
    mantid_config_path = os.path.join(dotmantid, 'Mantid.user.properties')
    with open(mantid_config_path, 'wt') as of:
        of.write('default.facility=%s\n' % facility)
        of.write('default.instrument=%s\n' % instrument)
        CORES = CORES or os.environ.get('CORES') or 2
        of.write('MultiThreaded.MaxCores = %s\n' % CORES)
    return
# this should be done before mantid is imported
_createDefaultMantidUserConfig()
