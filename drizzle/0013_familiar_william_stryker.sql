CREATE TABLE `cheminformatics_results` (
	`id` int AUTO_INCREMENT NOT NULL,
	`user_id` int NOT NULL,
	`input_smiles` text NOT NULL,
	`canonical_smiles` text,
	`workflow` enum('similarity','brics','validate','full_pipeline') NOT NULL,
	`threshold` decimal(3,2) DEFAULT '0.70',
	`max_hits` int DEFAULT 25,
	`results` json NOT NULL,
	`execution_time` int,
	`status` enum('pending','running','completed','failed') DEFAULT 'completed',
	`notes` text,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `cheminformatics_results_id` PRIMARY KEY(`id`)
);
--> statement-breakpoint
ALTER TABLE `cheminformatics_results` ADD CONSTRAINT `cheminformatics_results_user_id_users_id_fk` FOREIGN KEY (`user_id`) REFERENCES `users`(`id`) ON DELETE no action ON UPDATE no action;