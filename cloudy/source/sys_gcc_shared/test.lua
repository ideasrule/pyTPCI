-- Simple test case for use of Cloudy from the Lua scripting language
require 'cloudy'
local outfile='test_lua.out'
cloudy.init()
cloudy.output(outfile,'w') -- set Cloudy's output file
cloudy.read('test')
assert(cloudy.drive())
cloudy.output('','')              -- set output to stdout, close other stream
print('Finished test successfully -- see '..outfile..' for results')
