from fastapi import FastAPI, Depends, HTTPException
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from pydantic import BaseModel
from typing import Optional
from jose import JWTError, jwt
from passlib.context import CryptContext
from enum import Enum
import os
from datetime import datetime, timedelta

class UserRole(str, Enum):
    RESEARCHER = "researcher"
    ADMIN = "admin"
    VIEWER = "viewer"

class Token(BaseModel):
    access_token: str
    token_type: str

class TokenData(BaseModel):
    username: Optional[str] = None
    role: Optional[UserRole] = None

# In-memory user database for demonstration
FAKE_USERS_DB = {
    "researcher": {
        "username": "researcher",
        "full_name": "Researcher",
        "email": "researcher@example.com",
        "hashed_password": "$2b$12$EixZaB6Cq.w.t.h.i.s.I.s.A.Fake.Hash.1", # "password"
        "role": UserRole.RESEARCHER,
        "disabled": False,
    },
    "admin": {
        "username": "admin",
        "full_name": "Admin",
        "email": "admin@example.com",
        "hashed_password": "$2b$12$EixZaB6Cq.w.t.h.i.s.I.s.A.Fake.Hash.2", # "password"
        "role": UserRole.ADMIN,
        "disabled": False,
    }
}

class SecurityManager:
    def __init__(self):
        self.pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")
        self.secret_key = os.getenv("SECRET_KEY", "mysecretkey")
        self.algorithm = "HS256"
        self.oauth2_scheme = OAuth2PasswordBearer(tokenUrl="auth-service/token")

    def verify_password(self, plain_password, hashed_password):
        return self.pwd_context.verify(plain_password, hashed_password)

    def get_password_hash(self, password):
        return self.pwd_context.hash(password)

    def create_access_token(self, data: dict, expires_delta: timedelta = None):
        to_encode = data.copy()
        if expires_delta:
            expire = datetime.utcnow() + expires_delta
        else:
            expire = datetime.utcnow() + timedelta(minutes=15)
        to_encode.update({"exp": expire})
        encoded_jwt = jwt.encode(to_encode, self.secret_key, algorithm=self.algorithm)
        return encoded_jwt

    def get_current_user(self, token: str = Depends(lambda self: self.oauth2_scheme)):
        credentials_exception = HTTPException(
            status_code=401,
            detail="Could not validate credentials",
            headers={"WWW-Authenticate": "Bearer"},
        )
        try:
            payload = jwt.decode(token, self.secret_key, algorithms=[self.algorithm])
            username: str = payload.get("sub")
            role: UserRole = payload.get("role")
            if username is None:
                raise credentials_exception
            token_data = TokenData(username=username, role=role)
        except JWTError:
            raise credentials_exception
        
        user = FAKE_USERS_DB.get(token_data.username)
        if user is None or user["disabled"]:
            raise credentials_exception
        return user

    def verify_permission(self, required_role: UserRole, user_role: UserRole) -> bool:
        role_hierarchy = {
            UserRole.VIEWER: 1,
            UserRole.RESEARCHER: 2,
            UserRole.ADMIN: 3
        }
        return role_hierarchy.get(user_role, 0) >= role_hierarchy.get(required_role, 0)

app = FastAPI()
security_manager = SecurityManager()

@app.post("/token", response_model=Token)
async def login_for_access_token(form_data: OAuth2PasswordRequestForm = Depends()):
    user = FAKE_USERS_DB.get(form_data.username)
    if not user or not security_manager.verify_password(form_data.password, user["hashed_password"]):
        raise HTTPException(
            status_code=401,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token_expires = timedelta(minutes=30)
    access_token = security_manager.create_access_token(
        data={"sub": user["username"], "role": user["role"]}, expires_delta=access_token_expires
    )
    return {"access_token": access_token, "token_type": "bearer"}

@app.get("/users/me")
async def read_users_me(current_user: dict = Depends(lambda: security_manager.get_current_user(Depends(security_manager.oauth2_scheme)))):
    return current_user
