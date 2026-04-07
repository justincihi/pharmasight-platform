CREATE TABLE `metabolites` (
	`id` int AUTO_INCREMENT NOT NULL,
	`parent_analog_id` int NOT NULL,
	`smiles` text NOT NULL,
	`transformation` varchar(255) NOT NULL,
	`phase` enum('Phase I','Phase II') NOT NULL,
	`enzyme` varchar(64) NOT NULL,
	`probability` varchar(20) NOT NULL,
	`molecular_weight` varchar(20),
	`log_p` varchar(20),
	`metabolic_stability_score` int,
	`admet_score` int,
	`docking_score` int,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `metabolites_id` PRIMARY KEY(`id`)
);
