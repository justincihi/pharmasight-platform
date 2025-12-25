# PharmaSight Admin Dashboard

This directory contains the administrative dashboard for PharmaSight Platform.

## Architecture

The admin dashboard is a modern TypeScript/React application that provides:
- User management and authentication
- Database administration interface
- Analytics and reporting
- System configuration
- API key management for LLM integrations

## Technology Stack

- **Frontend**: React + Vite + TailwindCSS
- **Backend**: Node.js + Express
- **Database**: PostgreSQL with Drizzle ORM
- **Build**: TypeScript + pnpm

## Directory Structure

```
admin-dashboard/
├── client/          # React frontend application
├── server/          # Node.js backend API
├── drizzle/         # Database migrations
├── package.json     # Dependencies
└── vite.config.ts   # Build configuration
```

## Setup

### Prerequisites
- Node.js 18+
- pnpm
- PostgreSQL database

### Installation

```bash
cd admin-dashboard
pnpm install
```

### Configuration

Create a `.env` file in this directory:

```env
DATABASE_URL=postgresql://user:password@localhost:5432/pharmasight
FLASK_API_URL=http://localhost:5000
JWT_SECRET=your-secret-key
PORT=3000
```

### Development

```bash
# Start development server
pnpm dev

# Build for production
pnpm build

# Run production server
pnpm start
```

## Integration with Main Application

The admin dashboard connects to the main PharmaSight Flask application via REST API.

### Shared Database
Both applications use the same PostgreSQL database:
- Flask app uses SQLAlchemy ORM
- Admin dashboard uses Drizzle ORM
- Shared tables for users, compounds, research data

### API Integration
- Admin dashboard makes API calls to Flask backend
- JWT authentication for secure communication
- CORS configured for subdomain access

### Deployment Architecture

```
pharmasight.com
├── app.pharmasight.com     → Main Flask application (port 5000)
└── admin.pharmasight.com   → Admin dashboard (port 3000)
```

## Features

### User Management
- Create, read, update, delete users
- Role-based access control (RBAC)
- Session management

### Database Administration
- View and edit compound data
- Manage research findings
- Export data in various formats

### LLM API Integration
- Configure API keys (OpenAI, Perplexity, Gemini)
- Test API connections
- Monitor usage and costs

### Analytics
- Platform usage statistics
- Research activity tracking
- Performance metrics

## API Endpoints

The admin dashboard exposes its own REST API:

- `GET /api/users` - List all users
- `POST /api/users` - Create new user
- `GET /api/compounds` - List compounds
- `POST /api/config/llm` - Configure LLM API keys
- `GET /api/analytics` - Get platform analytics

## Security

- JWT-based authentication
- Environment variable for secrets
- HTTPS required in production
- CORS restricted to known origins
- Rate limiting on API endpoints

## Development Notes

- Hot module replacement (HMR) enabled in dev mode
- TypeScript strict mode enabled
- ESLint and Prettier configured
- Vitest for unit testing

## Production Deployment

### Using Docker

```bash
# Build Docker image
docker build -t pharmasight-admin .

# Run container
docker run -p 3000:3000 \
  -e DATABASE_URL=postgresql://... \
  -e FLASK_API_URL=http://app:5000 \
  pharmasight-admin
```

### Using PM2

```bash
# Build application
pnpm build

# Start with PM2
pm2 start npm --name "pharmasight-admin" -- start
```

## Troubleshooting

### Database Connection Issues
- Verify DATABASE_URL is correct
- Check PostgreSQL is running
- Ensure database exists and migrations are applied

### API Connection Issues
- Verify FLASK_API_URL points to Flask app
- Check CORS configuration
- Ensure both apps can communicate (network/firewall)

### Build Issues
- Clear node_modules and reinstall: `rm -rf node_modules && pnpm install`
- Clear build cache: `rm -rf dist .vite`
- Check Node.js version: `node --version` (should be 18+)

## Contributing

When making changes to the admin dashboard:
1. Create a feature branch
2. Make changes and test locally
3. Run tests: `pnpm test`
4. Build to verify: `pnpm build`
5. Commit and push changes

## License

Proprietary - PharmaSight Platform
