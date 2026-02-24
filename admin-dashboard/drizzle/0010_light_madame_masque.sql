ALTER TABLE `analog_discoveries` ADD `parent_analog_id` int;--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `optimization_generation` int DEFAULT 1;--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `optimization_target` varchar(128);--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `optimization_notes` text;--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `toxicity_profile` text;--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `synthetic_accessibility` text;--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `metabolites` text;