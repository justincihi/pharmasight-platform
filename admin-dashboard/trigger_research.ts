import { runSchedulerNow } from './server/autonomousScheduler';

console.log('🚀 Starting manual autonomous research run...');
console.log('Time:', new Date().toISOString());

runSchedulerNow()
  .then(() => {
    console.log('✅ Manual run completed successfully');
    process.exit(0);
  })
  .catch((err) => {
    console.error('❌ Manual run failed:', err);
    process.exit(1);
  });
