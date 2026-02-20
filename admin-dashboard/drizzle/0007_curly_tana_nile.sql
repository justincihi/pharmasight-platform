CREATE TABLE `docking_queue` (
	`id` int AUTO_INCREMENT NOT NULL,
	`analog_id` int NOT NULL,
	`target` varchar(100) NOT NULL,
	`status` enum('pending','running','completed','failed') NOT NULL DEFAULT 'pending',
	`priority` int NOT NULL DEFAULT 5,
	`binding_affinity` varchar(50),
	`docking_score` int,
	`ligand_pdb` text,
	`receptor_pdb` text,
	`error_message` text,
	`started_at` timestamp,
	`completed_at` timestamp,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `docking_queue_id` PRIMARY KEY(`id`)
);
