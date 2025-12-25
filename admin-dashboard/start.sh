#!/bin/bash
# Start script for Admin Dashboard on Render

set -e

echo "Starting PharmaSight Admin Dashboard..."

# Install dependencies if not already installed
if [ ! -d "node_modules" ]; then
    echo "Installing dependencies..."
    pnpm install --frozen-lockfile
fi

# Build the application
echo "Building application..."
pnpm build

# Start the production server
echo "Starting server..."
exec pnpm start
