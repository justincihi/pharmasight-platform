CREATE TABLE `batch_docking_jobs` (
	`id` varchar(64) NOT NULL,
	`user_id` varchar(128) NOT NULL,
	`job_name` varchar(255) NOT NULL,
	`status` enum('pending','running','completed','failed','cancelled') NOT NULL DEFAULT 'pending',
	`total_compounds` int NOT NULL,
	`completed_compounds` int NOT NULL DEFAULT 0,
	`failed_compounds` int NOT NULL DEFAULT 0,
	`target_name` varchar(255) NOT NULL,
	`parameters_id` int,
	`results_summary` text,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`started_at` timestamp,
	`completed_at` timestamp,
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `batch_docking_jobs_id` PRIMARY KEY(`id`)
);
--> statement-breakpoint
CREATE TABLE `batch_docking_results` (
	`id` int AUTO_INCREMENT NOT NULL,
	`job_id` varchar(64) NOT NULL,
	`analog_id` int NOT NULL,
	`status` enum('pending','completed','failed') NOT NULL DEFAULT 'pending',
	`binding_affinity` varchar(64),
	`docking_score` int,
	`num_poses` int,
	`top_poses` text,
	`error_message` text,
	`completed_at` timestamp,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `batch_docking_results_id` PRIMARY KEY(`id`)
);
--> statement-breakpoint
CREATE TABLE `docking_parameters` (
	`id` int AUTO_INCREMENT NOT NULL,
	`user_id` varchar(128) NOT NULL,
	`name` varchar(255) NOT NULL,
	`target_name` varchar(255) NOT NULL,
	`box_center_x` varchar(64) NOT NULL,
	`box_center_y` varchar(64) NOT NULL,
	`box_center_z` varchar(64) NOT NULL,
	`box_size_x` varchar(64) NOT NULL,
	`box_size_y` varchar(64) NOT NULL,
	`box_size_z` varchar(64) NOT NULL,
	`exhaustiveness` int NOT NULL DEFAULT 8,
	`num_poses` int NOT NULL DEFAULT 9,
	`is_default` tinyint DEFAULT 0,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `docking_parameters_id` PRIMARY KEY(`id`)
);
--> statement-breakpoint
CREATE TABLE `pdb_receptors` (
	`id` varchar(64) NOT NULL,
	`name` varchar(255) NOT NULL,
	`file_key` varchar(512) NOT NULL,
	`url` text NOT NULL,
	`uploaded_by` varchar(128) NOT NULL,
	`uploaded_at` timestamp NOT NULL DEFAULT (now()),
	`file_size` int NOT NULL,
	`target_name` varchar(255) NOT NULL,
	`description` text,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `pdb_receptors_id` PRIMARY KEY(`id`)
);
