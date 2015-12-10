#!/usr/bin/env th

require 'hdf5'

require 'convnet'

----------------------------------------------------------------
-- parse arguments
----------------------------------------------------------------
cmd = torch.CmdLine()
cmd:text()
cmd:text('DNA ConvNet Second Layer Weights')
cmd:text()
cmd:text('Arguments')
cmd:argument('model_file')
cmd:argument('out_file')
cmd:text()
opt = cmd:parse(arg)

-- fix seed
torch.manualSeed(1)

----------------------------------------------------------------
-- load data
----------------------------------------------------------------
local convnet_params = torch.load(opt.model_file)
-- print('fml')
-- for k,v in pairs(convnet_params) do
-- 	print(v)
-- en

local convnet = ConvNet:__init()
convnet:load(convnet_params)
print(convnet)

print ("modules")
for k,v in pairs(convnet.model.modules) do
	print (k,v)
end

-- get convolution filter params
local filter_weights = convnet.model.modules[5].weight:squeeze()

-- dump to file, load into python.
local hdf_out = hdf5.open(opt.out_file, 'w')
hdf_out:write('weights', filter_weights)
hdf_out:close()