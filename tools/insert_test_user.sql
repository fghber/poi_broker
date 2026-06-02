-- Insert a validated test user into users.db for testing and local development purposes ONLY! 
-- DO NOT USE THIS ACCOUNT FOR ANY REAL PURPOSES, AND DO NOT SHARE THE CREDENTIALS WITH ANYONE OUTSIDE OF A SAFE TESTING ENVIRONMENT.
-- Password: testpass123  (werkzeug scrypt hash)
-- email_verified = 1 bypasses email confirmation requirement
--
-- Note: email must pass RFC 5321 syntax (use a real-looking domain such as example.com)
--
-- Usage:
--   sqlite3 /path/to/users.db < tools/insert_test_user.sql
--
-- To remove the test user:
--   DELETE FROM user WHERE email = 'testuser@example.com';

INSERT OR IGNORE INTO user (
    email,
    password,
    name,
    role,
    reset_token,
    reset_token_expires,
    email_verified,
    email_verification_token
) VALUES (
    'testuser@example.com',
    'scrypt:32768:8:1$VHh06SV6TGKitB3X$936f8348fcbbc52d2294ea1ffc07107341a101013aef2697343c4be4043c15db213c43807b9f8fab0b85cf8529f857f55dda19345d9df1a968f2be9126971cd7',
    'Test User',
    'user',
    NULL,
    NULL,
    1,
    NULL
);
