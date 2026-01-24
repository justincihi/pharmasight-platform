CREATE TABLE `analog_discoveries` (
	`id` int AUTO_INCREMENT NOT NULL,
	`compound_id` varchar(128) NOT NULL,
	`compound_name` varchar(255) NOT NULL,
	`parent_compound` varchar(255) NOT NULL,
	`smiles` text NOT NULL,
	`confidence_score` int NOT NULL,
	`similarity_score` int NOT NULL,
	`safety_score` int NOT NULL,
	`efficacy_score` int NOT NULL,
	`drug_likeness_score` int NOT NULL,
	`patent_status` enum('patent-free','patent-opportunity','patented','unknown') NOT NULL,
	`patent_numbers` text,
	`fda_status` varchar(64),
	`market_value` varchar(64),
	`therapeutic_potential` text,
	`key_differences` text,
	`mechanism_of_action` text,
	`molecular_weight` varchar(64),
	`log_p` varchar(64),
	`h_bond_donors` int,
	`h_bond_acceptors` int,
	`pubchem_cid` varchar(64),
	`chembl_id` varchar(64),
	`discovered_by` varchar(128) NOT NULL,
	`discovery_method` varchar(128),
	`discovered_at` timestamp NOT NULL DEFAULT (now()),
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `analog_discoveries_id` PRIMARY KEY(`id`),
	CONSTRAINT `analog_discoveries_compound_id_unique` UNIQUE(`compound_id`)
);
--> statement-breakpoint
CREATE TABLE `notifications` (
	`id` int AUTO_INCREMENT NOT NULL,
	`user_id` int NOT NULL,
	`analog_id` int,
	`title` varchar(255) NOT NULL,
	`message` text NOT NULL,
	`notification_type` enum('new-discovery','high-confidence','patent-alert','system') NOT NULL,
	`is_read` int NOT NULL DEFAULT 0,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `notifications_id` PRIMARY KEY(`id`)
);
--> statement-breakpoint
CREATE TABLE `test_results` (
	`id` int AUTO_INCREMENT NOT NULL,
	`analog_id` int NOT NULL,
	`test_type` enum('admet','docking','toxicity','pkpd','quantum') NOT NULL,
	`test_status` enum('pending','running','completed','failed') NOT NULL,
	`results` text,
	`error_message` text,
	`run_by` int NOT NULL,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`completed_at` timestamp,
	CONSTRAINT `test_results_id` PRIMARY KEY(`id`)
);
