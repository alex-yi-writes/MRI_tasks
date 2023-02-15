function [tmplist] = makeRetrievalList(path_parent,ID)
%% set paths

path_model = [path_parent 'stim/stiminfo/models/'];
load([path_model eval(['REPLAY_RET7T_GA_designstruct_block' num2str(BlockNum) '.mat'])])





end