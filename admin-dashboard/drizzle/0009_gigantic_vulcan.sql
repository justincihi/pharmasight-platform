CREATE TABLE `bookmarks` (
	`id` int AUTO_INCREMENT NOT NULL,
	`user_id` int NOT NULL,
	`analog_id` int,
	`notification_id` int,
	`title` varchar(255) NOT NULL,
	`notes` text,
	`category` enum('high-priority','review-later','promising','archived') NOT NULL DEFAULT 'review-later',
	`created_at` timestamp NOT NULL DEFAULT (now()),
	`updated_at` timestamp NOT NULL DEFAULT (now()) ON UPDATE CURRENT_TIMESTAMP,
	CONSTRAINT `bookmarks_id` PRIMARY KEY(`id`)
);
