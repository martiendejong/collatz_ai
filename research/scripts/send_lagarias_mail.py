"""One-shot: send the Lagarias record mail from info@martiendejong.nl.
Credential: vault.prospergenics.com project 8, credential 20 (IMAP/SMTP).
Body: taken verbatim from lagarias-mail-draft.md (English part only).
Bcc to info@martiendejong.nl so Martien has the sent copy in his inbox.
"""
import requests, smtplib, ssl, sys, re
from email.mime.text import MIMEText
from email.utils import formatdate, make_msgid

VAULT_KEY = open(r"C:\home\claude\.claude\vault\vault.key").read().strip()
h = {"X-API-Key": VAULT_KEY}
cred = requests.get(
    "https://vault.prospergenics.com/api/projects/8/credentials/20",
    headers=h, timeout=30).json()
user = cred.get("username") or "info@martiendejong.nl"
pw = cred["password"]
print(f"credential ok (user: {user})")

draft = open(r"E:\projects\collatz\research\lagarias-mail-draft.md",
             encoding="utf-8").read()
m = re.search(r"\*\*Subject: (.+?)\*\*\n\n(Dear Professor.*?info@martiendejong\.nl)\n",
              draft, re.S)
subject = m.group(1)
body = m.group(2)
print(f"subject: {subject}")
print(f"body: {len(body)} chars, starts: {body[:40]!r}")

TO = "lagarias@umich.edu"
FROM = "info@martiendejong.nl"
BCC = "info@martiendejong.nl"

msg = MIMEText(body, "plain", "utf-8")
msg["Subject"] = subject
msg["From"] = f"Martien de Jong <{FROM}>"
msg["To"] = TO
msg["Date"] = formatdate(localtime=True)
msg["Message-ID"] = make_msgid(domain="martiendejong.nl")

if "--send" not in sys.argv:
    print("\nDRY RUN (pass --send to transmit). Full message below:\n")
    print(msg.as_string()[:600])
    sys.exit(0)

ctx = ssl.create_default_context()
try:
    s = smtplib.SMTP_SSL("mail.martiendejong.nl", 465, context=ctx, timeout=60)
except Exception:
    s = smtplib.SMTP("mail.martiendejong.nl", 587, timeout=60)
    s.starttls(context=ctx)
s.login(user, pw)
s.sendmail(FROM, [TO, BCC], msg.as_string())
s.quit()
print(f"SENT to {TO} (bcc {BCC})")
