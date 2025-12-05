#!/usr/bin/env python3
"""
PharmaSight™ Database-Backed Authentication Module
Secure user authentication with hashed passwords stored in PostgreSQL
Includes TOTP-based Two-Factor Authentication (2FA)
"""

import os
import psycopg2
from psycopg2.extras import RealDictCursor
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin
from datetime import datetime
import secrets
import pyotp
import base64
from io import BytesIO

DATABASE_URL = os.environ.get('DATABASE_URL')

class User(UserMixin):
    """User model for Flask-Login"""
    def __init__(self, id, username, email, role, created_at, last_login=None):
        self.id = id
        self.username = username
        self.email = email
        self.role = role
        self.created_at = created_at
        self.last_login = last_login
    
    def get_id(self):
        return str(self.id)
    
    @property
    def is_admin(self):
        return self.role == 'admin'


def get_db_connection():
    """Get a database connection"""
    return psycopg2.connect(DATABASE_URL, cursor_factory=RealDictCursor)


def init_auth_tables():
    """Initialize authentication tables in the database"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    cur.execute('''
        CREATE TABLE IF NOT EXISTS users (
            id SERIAL PRIMARY KEY,
            username VARCHAR(50) UNIQUE NOT NULL,
            email VARCHAR(255) UNIQUE NOT NULL,
            password_hash VARCHAR(255) NOT NULL,
            role VARCHAR(20) DEFAULT 'user',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            last_login TIMESTAMP,
            is_active BOOLEAN DEFAULT TRUE
        )
    ''')
    
    cur.execute('''
        CREATE TABLE IF NOT EXISTS login_attempts (
            id SERIAL PRIMARY KEY,
            username VARCHAR(50),
            ip_address VARCHAR(45),
            success BOOLEAN,
            attempted_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    cur.execute('''
        CREATE TABLE IF NOT EXISTS user_sessions (
            id SERIAL PRIMARY KEY,
            user_id INTEGER REFERENCES users(id),
            session_token VARCHAR(255) UNIQUE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            expires_at TIMESTAMP,
            is_valid BOOLEAN DEFAULT TRUE
        )
    ''')
    
    cur.execute('''
        CREATE TABLE IF NOT EXISTS user_totp (
            id SERIAL PRIMARY KEY,
            user_id INTEGER REFERENCES users(id) UNIQUE,
            totp_secret VARCHAR(64) NOT NULL,
            is_enabled BOOLEAN DEFAULT FALSE,
            backup_codes TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            verified_at TIMESTAMP
        )
    ''')
    
    conn.commit()
    cur.close()
    conn.close()
    print("Authentication tables initialized successfully")


def create_user(username, email, password, role='user'):
    """Create a new user with hashed password"""
    password_hash = generate_password_hash(password, method='pbkdf2:sha256')
    
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('''
            INSERT INTO users (username, email, password_hash, role)
            VALUES (%s, %s, %s, %s)
            RETURNING id, username, email, role, created_at
        ''', (username, email, password_hash, role))
        
        user_data = cur.fetchone()
        conn.commit()
        
        return {
            'success': True,
            'user': {
                'id': user_data['id'],
                'username': user_data['username'],
                'email': user_data['email'],
                'role': user_data['role']
            }
        }
    except psycopg2.IntegrityError as e:
        conn.rollback()
        if 'username' in str(e):
            return {'success': False, 'error': 'Username already exists'}
        elif 'email' in str(e):
            return {'success': False, 'error': 'Email already exists'}
        return {'success': False, 'error': 'User creation failed'}
    finally:
        cur.close()
        conn.close()


def authenticate_user(username, password, ip_address=None):
    """Authenticate a user and return user object if successful"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('''
            SELECT id, username, email, password_hash, role, created_at, last_login, is_active
            FROM users
            WHERE username = %s OR email = %s
        ''', (username, username))
        
        user_data = cur.fetchone()
        
        if user_data and user_data['is_active'] and check_password_hash(user_data['password_hash'], password):
            cur.execute('''
                UPDATE users SET last_login = CURRENT_TIMESTAMP WHERE id = %s
            ''', (user_data['id'],))
            
            cur.execute('''
                INSERT INTO login_attempts (username, ip_address, success)
                VALUES (%s, %s, %s)
            ''', (username, ip_address, True))
            
            conn.commit()
            
            return User(
                id=user_data['id'],
                username=user_data['username'],
                email=user_data['email'],
                role=user_data['role'],
                created_at=user_data['created_at'],
                last_login=datetime.now()
            )
        else:
            cur.execute('''
                INSERT INTO login_attempts (username, ip_address, success)
                VALUES (%s, %s, %s)
            ''', (username, ip_address, False))
            conn.commit()
            return None
            
    finally:
        cur.close()
        conn.close()


def get_user_by_id(user_id):
    """Get user by ID for Flask-Login"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('''
            SELECT id, username, email, role, created_at, last_login
            FROM users WHERE id = %s AND is_active = TRUE
        ''', (user_id,))
        
        user_data = cur.fetchone()
        
        if user_data:
            return User(
                id=user_data['id'],
                username=user_data['username'],
                email=user_data['email'],
                role=user_data['role'],
                created_at=user_data['created_at'],
                last_login=user_data['last_login']
            )
        return None
    finally:
        cur.close()
        conn.close()


def change_password(user_id, old_password, new_password):
    """Change user password"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('SELECT password_hash FROM users WHERE id = %s', (user_id,))
        user_data = cur.fetchone()
        
        if not user_data or not check_password_hash(user_data['password_hash'], old_password):
            return {'success': False, 'error': 'Current password is incorrect'}
        
        new_hash = generate_password_hash(new_password, method='pbkdf2:sha256')
        cur.execute('UPDATE users SET password_hash = %s WHERE id = %s', (new_hash, user_id))
        conn.commit()
        
        return {'success': True, 'message': 'Password changed successfully'}
    finally:
        cur.close()
        conn.close()


def get_all_users():
    """Get all users (admin only)"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('''
            SELECT id, username, email, role, created_at, last_login, is_active
            FROM users ORDER BY created_at DESC
        ''')
        
        users = cur.fetchall()
        return [dict(u) for u in users]
    finally:
        cur.close()
        conn.close()


def delete_user(user_id):
    """Soft delete a user"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('UPDATE users SET is_active = FALSE WHERE id = %s', (user_id,))
        conn.commit()
        return {'success': True}
    finally:
        cur.close()
        conn.close()


def create_default_admin():
    """Create default admin user if no users exist"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('SELECT COUNT(*) as count FROM users')
        count = cur.fetchone()['count']
        
        if count == 0:
            admin_password = os.environ.get('ADMIN_PASSWORD', secrets.token_urlsafe(16))
            result = create_user(
                username='admin',
                email='admin@pharmasight.local',
                password=admin_password,
                role='admin'
            )
            
            if result['success']:
                print(f"\n{'='*60}")
                print("DEFAULT ADMIN ACCOUNT CREATED")
                print(f"{'='*60}")
                print(f"Username: admin")
                print(f"Password: {admin_password}")
                print(f"{'='*60}")
                print("IMPORTANT: Change this password immediately!")
                print(f"{'='*60}\n")
                return admin_password
        return None
    finally:
        cur.close()
        conn.close()


def reset_user_password(user_id, new_password=None):
    """Admin function to reset a user's password"""
    if new_password is None:
        new_password = secrets.token_urlsafe(12)
    
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        new_hash = generate_password_hash(new_password, method='pbkdf2:sha256')
        cur.execute('UPDATE users SET password_hash = %s WHERE id = %s', (new_hash, user_id))
        conn.commit()
        
        return {'success': True, 'new_password': new_password}
    finally:
        cur.close()
        conn.close()


# ========== TOTP 2FA FUNCTIONS ==========

def generate_totp_secret(user_id):
    """Generate a new TOTP secret for a user and store it"""
    totp_secret = pyotp.random_base32()
    
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        # Generate backup codes
        backup_codes = [secrets.token_hex(4).upper() for _ in range(8)]
        backup_codes_str = ','.join(backup_codes)
        
        # Insert or update TOTP secret
        cur.execute('''
            INSERT INTO user_totp (user_id, totp_secret, backup_codes, is_enabled)
            VALUES (%s, %s, %s, FALSE)
            ON CONFLICT (user_id) 
            DO UPDATE SET totp_secret = %s, backup_codes = %s, is_enabled = FALSE, verified_at = NULL
        ''', (user_id, totp_secret, backup_codes_str, totp_secret, backup_codes_str))
        
        conn.commit()
        
        return {
            'success': True,
            'secret': totp_secret,
            'backup_codes': backup_codes
        }
    except Exception as e:
        conn.rollback()
        return {'success': False, 'error': str(e)}
    finally:
        cur.close()
        conn.close()


def get_totp_provisioning_uri(user_id, username):
    """Get the provisioning URI for TOTP setup (for QR code generation)"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('SELECT totp_secret FROM user_totp WHERE user_id = %s', (user_id,))
        result = cur.fetchone()
        
        if not result:
            return None
        
        totp = pyotp.TOTP(result['totp_secret'])
        uri = totp.provisioning_uri(name=username, issuer_name='PharmaSight')
        
        return uri
    finally:
        cur.close()
        conn.close()


def generate_totp_qr_code(user_id, username):
    """Generate a QR code image for TOTP setup"""
    try:
        import qrcode
        
        uri = get_totp_provisioning_uri(user_id, username)
        if not uri:
            return None
        
        # Generate QR code
        qr = qrcode.QRCode(version=1, box_size=10, border=5)
        qr.add_data(uri)
        qr.make(fit=True)
        
        img = qr.make_image(fill_color="black", back_color="white")
        
        # Convert to base64
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_base64 = base64.b64encode(buffered.getvalue()).decode()
        
        return f"data:image/png;base64,{img_base64}"
    except Exception as e:
        print(f"Error generating QR code: {e}")
        return None


def verify_totp_code(user_id, code):
    """Verify a TOTP code for a user"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('SELECT totp_secret, backup_codes FROM user_totp WHERE user_id = %s', (user_id,))
        result = cur.fetchone()
        
        if not result:
            return {'success': False, 'error': 'TOTP not configured'}
        
        totp_secret = result['totp_secret']
        backup_codes = result['backup_codes'].split(',') if result['backup_codes'] else []
        
        # Check TOTP code
        totp = pyotp.TOTP(totp_secret)
        if totp.verify(code, valid_window=1):  # Allow 30 seconds window
            return {'success': True, 'method': 'totp'}
        
        # Check backup codes
        if code.upper() in backup_codes:
            # Remove used backup code
            backup_codes.remove(code.upper())
            cur.execute(
                'UPDATE user_totp SET backup_codes = %s WHERE user_id = %s',
                (','.join(backup_codes), user_id)
            )
            conn.commit()
            return {'success': True, 'method': 'backup_code', 'remaining_codes': len(backup_codes)}
        
        return {'success': False, 'error': 'Invalid code'}
    finally:
        cur.close()
        conn.close()


def enable_totp(user_id, code):
    """Enable TOTP after verifying the initial code"""
    verify_result = verify_totp_code(user_id, code)
    
    if not verify_result['success']:
        return verify_result
    
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('''
            UPDATE user_totp 
            SET is_enabled = TRUE, verified_at = CURRENT_TIMESTAMP 
            WHERE user_id = %s
        ''', (user_id,))
        conn.commit()
        
        return {'success': True, 'message': 'Two-factor authentication enabled'}
    finally:
        cur.close()
        conn.close()


def disable_totp(user_id):
    """Disable TOTP for a user"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('DELETE FROM user_totp WHERE user_id = %s', (user_id,))
        conn.commit()
        return {'success': True, 'message': 'Two-factor authentication disabled'}
    finally:
        cur.close()
        conn.close()


def get_totp_status(user_id):
    """Check if TOTP is enabled for a user"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        cur.execute('''
            SELECT is_enabled, verified_at, 
                   CASE WHEN backup_codes IS NOT NULL THEN 
                       array_length(string_to_array(backup_codes, ','), 1) 
                   ELSE 0 END as backup_codes_remaining
            FROM user_totp 
            WHERE user_id = %s
        ''', (user_id,))
        result = cur.fetchone()
        
        if not result:
            return {'enabled': False, 'configured': False}
        
        return {
            'enabled': result['is_enabled'],
            'configured': True,
            'verified_at': result['verified_at'].isoformat() if result['verified_at'] else None,
            'backup_codes_remaining': result['backup_codes_remaining'] or 0
        }
    finally:
        cur.close()
        conn.close()


def is_totp_required(user_id):
    """Check if TOTP verification is required for login"""
    status = get_totp_status(user_id)
    return status.get('enabled', False)
