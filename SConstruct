
env = Environment(CCFLAGS='-Wall') # show all warnings

program = env.Program(target='csufmedia_test',
                      source=['csufmedia.cpp', 'csufmedia_test.cpp'],
                      LIBS=['boost_unit_test_framework'],
                      LINKFLAGS='--static') # must link statically to workaround boost bug

# Always run unit tests, copied from http://www.scons.org/wiki/UnitTests
# must be env.Alias as discussed here:
# http://scons.tigris.org/ds/viewMessage.do?dsForumId=1271&dsMessageId=2459175&orderBy=createDate&orderType=desc
test_alias = env.Alias('test', [program], program[0].abspath)
AlwaysBuild(test_alias)
