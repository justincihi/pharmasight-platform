ALTER TABLE `analog_discoveries` ADD `approval_status` enum('pending','approved','rejected') DEFAULT 'pending' NOT NULL;--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `approved_by` varchar(128);--> statement-breakpoint
ALTER TABLE `analog_discoveries` ADD `approved_at` timestamp;