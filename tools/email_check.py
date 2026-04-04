import smtplib
import ssl
import socket
import sys
 
# --- CONFIGURATION ---
SMTP_HOST = 'smtp.gmail.com'
SMTP_PORT = 587 # Using 587 for STARTTLS as it's better for debugging
USER = 'your-new-email@gmail.com'
PASS = 'xxxx xxxx xxxx xxxx' # Your 16-digit App Password (no spaces)
DEST = 'your-personal-email@example.com'
 
def run_debug_test():
    print(f"--- Starting SMTP Debug Test (Target: {SMTP_HOST}) ---")
 
    try:
        # 1. Force IPv4 Resolution
        # print(f"[*] Resolving {SMTP_HOST} to IPv4...")
        # addr_info = socket.getaddrinfo(SMTP_HOST, SMTP_PORT, socket.AF_INET)
        # ipv4_address = addr_info[0][4][0]
        # print(f"[+] Using IPv4 Address: {ipv4_address}")
 
        # 2. Connect
        print("[*] Connecting to server...")
        server = smtplib.SMTP(SMTP_HOST, SMTP_PORT, timeout=10)
 
        # 3. Enable Full Debug Output
        # 0 = off, 1 = variable messages, 2 = all messages with timestamps
        server.set_debuglevel(1)
 
        # 4. Say Hello (EHLO)
        server.ehlo()
 
        # 5. Start TLS
        if server.has_extn('STARTTLS'):
            print("[*] Starting TLS encryption...")
            context = ssl.create_default_context()
            server.starttls(context=context)
            server.ehlo() # Re-identify after encryption
 
        # 6. Login
        print(f"[*] Attempting login for {USER}...")
        server.login(USER, PASS)
        print("[+] Login Successful!")
 
        # 7. Send a test ping
        msg = f"Subject: Server Test\n\nThis is a test from the Ubuntu server."
        server.sendmail(USER, DEST, msg)
        print("[+] Test email sent successfully!")
 
        server.quit()
 
    except socket.timeout:
        print("\n[!] ERROR: Connection timed out. Is Port 587 blocked by your provider (AWS/GCP/DigitalOcean)?")
    except smtplib.SMTPAuthenticationError as e:
        print(f"\n[!] AUTHENTICATION FAILED: {e.smtp_code} {e.smtp_error.decode()}")
        print(" Check: 1. App Password typo 2. DisplayUnlockCaptcha 3. Account Security Alerts")
    except Exception as e:
        print(f"\n[!] AN UNEXPECTED ERROR OCCURRED: {type(e).__name__}: {e}")
 
if __name__ == "__main__":
    if 'xxxx' in PASS:
        print("Error: Please update the USER and PASS variables in the script first.")
        sys.exit(1)
    run_debug_test()