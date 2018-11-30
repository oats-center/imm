clear all;

gd_events = imm_ule_identify(2018);
save('../mat-files/2018/combine_gd_ul.mat', 'gd_events');

clear all;

gd_events = imm_ule_identify(2017);
save('../mat-files/2017/combine_gd_ul.mat', 'gd_events');

clear all;

gd_events = imm_ule_identify(2016);
save('../mat-files/2016/combine_gd_ul.mat', 'gd_events');

clear all;

gd_events = imm_ule_identify(2015);
save('../mat-files/2015/combine_gd_ul.mat', 'gd_events');
