CREATE TABLE `analysis_results` (
	`id` int AUTO_INCREMENT NOT NULL,
	`analog_id` int NOT NULL,
	`analysis_type` enum('docking','toxicity','admet','pkpd') NOT NULL,
	`smiles` text NOT NULL,
	`target` varchar(128),
	`result` json NOT NULL,
	`source` enum('python','api','fallback') NOT NULL DEFAULT 'python',
	`execution_time` int,
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`created_by` varchar(64),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `analysis_results_id` PRIMARY KEY(`id`)
);
