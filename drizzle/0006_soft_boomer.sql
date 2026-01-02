ALTER TABLE `analog_discoveries` ADD `binding_affinity` varchar(64);--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `docking_score` int;--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `docking_target` varchar(128);