#!/usr/bin/env python3
# Script to test access to authenticated resources via REST interface.
# Written by Keith Jolley
# Copyright (c) 2017-2024, University of Oxford
# E-mail: keith.jolley@biology.ox.ac.uk
#
# This file is part of Bacterial Isolate Genome Sequence Database (BIGSdb).
#
# BIGSdb is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BIGSdb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# The test databases can be reached at https://pubmlst.org/test/.
# To use these, sign up for a PubMLST account (https://pubmlst.org/site_accounts.shtml)
# and link this account with the pubmlst_test_seqdef and pubmlst_test_isolates
# databases (https://pubmlst.org/site_accounts.shtml#registering_with_databases)
#
# Please note that the consumer key below will only work for access to the
# PubMLST test databases. Please contact keith.jolley@biology.ox.ac.uk to obtain your
# own consumer key for use in your own projects.
# Version 20241031

###
import argparse
import re
import os
import sys
import json
import base64
from rauth import OAuth1Service, OAuth1Session
from urllib.parse import parse_qs

from dotenv import load_dotenv
load_dotenv('/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/config_files/.env')

CONSUMER_KEY = os.getenv('CONSUMER_KEY')
CONSUMER_SECRET = os.getenv('CONSUMER_SECRET')
TEST_REST_URL = os.getenv('TEST_REST_URL')
SESSION_REST_URL = os.getenv('SESSION_REST_URL')
TEST_WEB_URL = os.getenv('TEST_WEB_URL')
ACCESS_TOKEN_PATH = os.getenv('ACCESS_TOKEN_PATH')
SESSION_TOKEN = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/config_files/.session_token'

REQUEST_TOKEN_URL = TEST_REST_URL + "/oauth/get_request_token"
ACCESS_TOKEN_URL = TEST_REST_URL + "/oauth/get_access_token"
AUTHORIZE_URL = TEST_WEB_URL + "&page=authorizeClient"

METHOD='POST'
ROUTE='/sequence'

# Access and session tokens are stored within current directory.
# If a session token is missing or expired, a new one will be requested using the access token.
# If an access token is missing or expired, a new one will be requested.
def main():

    service = OAuth1Service(
        name="BIGSdb",
        consumer_key=CONSUMER_KEY,
        consumer_secret=CONSUMER_SECRET,
        request_token_url=REQUEST_TOKEN_URL,
        access_token_url=ACCESS_TOKEN_URL,
        base_url=TEST_REST_URL,
    )

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--arguments",
        help="Data to put in during POST or additional GET arguments e.g. 'type=alleles&software=Enterobase'",
    )
    parser.add_argument("-f", "--file", help="Name of file to upload")
    parser.add_argument(
        "-i",
        "--isolates_file",
        help="Relative path of tab-delimited file of isolate data to upload",
    )
    parser.add_argument(
        "-reqm",
        "--request_method",
        help="Set HTTP method (default GET)",
        choices=["GET", "POST", "PUT", "DELETE"],
        default=METHOD,
    )
    parser.add_argument(
        "-p",
        "--profiles_file",
        help="Relative path of tab-delimited file of allelic profiles to upload",
    )
    parser.add_argument(
        "-r", "--route", help="Relative path of route, e.g. 'submissions'", default=ROUTE
    )
    parser.add_argument(
        "-s",
        "--sequence_file",
        help="Relative path of FASTA or single sequence file to upload",
        default=''
    )
    # parser.add_argument(
    #     "--prompt",
    #     help="Prompt before connection requests (used for demonstration purposes)",
    #     type=bool,
    # )
    args = parser.parse_args()

    (token, secret) = retrieve_token(SESSION_TOKEN)
    if not token or not secret:
        (token, secret) = get_session_token(None, None)
    get_route(args.route, token, secret, args)
    return


def retrieve_token(token_name):
    if not os.path.exists(token_name):
        return (None, None)
    file = open(token_name, "r")
    content = file.read()
    file.close()
    lines = content.split("\n")  # split it into lines
    token = ""
    secret = ""
    for line in lines:
        p = line.split("=")
        if p[0] == "token":
            token = p[1]
        if p[0] == "secret":
            secret = p[1]
    return (token, secret)


def trim_url_args(url):
    if not "?" in url:
        return url, {}
    trimmed_url, param_string = url.split("?")
    params = parse_qs(param_string)
    params = {k: int(v[0]) for k, v in params.items()}
    return trimmed_url, params


def get_route(route, token, secret, args):
    if route and not re.match(r"^/", route):
        route = "/" + route
    url = TEST_REST_URL + route
    print("Accessing authenticated resource ({0})...\n".format(url))
    session = OAuth1Session(
        CONSUMER_KEY, CONSUMER_SECRET, access_token=token, access_token_secret=secret
    )
    extra_params = {}
    if args.sequence_file:
        if not os.path.exists(args.sequence_file):
            print("Sequence file " + args.sequence_file + " does not exist.")
            sys.exit(1)
        with open(args.sequence_file, "r") as seq_file:
            data = seq_file.read()
            extra_params["sequence"] = data
    if args.profiles_file:
        if not os.path.exists(args.profiles_file):
            print("Profiles file " + args.profiles_file + " does not exist.")
            sys.exit(1)
        with open(args.profiles_file, "r") as profiles_file:
            data = profiles_file.read()
            extra_params["profiles"] = data
    if args.isolates_file:
        if not os.path.exists(args.isolates_file):
            print("Isolates file " + args.isolates_file + " does not exist.")
            sys.exit(1)
        with open(args.isolates_file, "r") as isolates_file:
            data = isolates_file.read()
            extra_params["isolates"] = data
    if args.file:
        if not os.path.exists(args.file):
            print("File " + args.file + " does not exist.")
            sys.exit(1)
        with open(args.file, mode="rb") as file:
            content = file.read()
            extra_params["upload"] = base64.b64encode(content)
    if args.arguments:
        p = args.arguments.split("&")
        for pa in p:
            ps = pa.split("=")
            extra_params[ps[0]] = ps[1]
    trimmed_url, request_params = trim_url_args(url)
    if args.request_method == "GET":
        r = session.get(trimmed_url, params=request_params)
    elif args.request_method == "POST":
        r = session.post(
            trimmed_url,
            params=request_params,
            data=json.dumps(extra_params),
            headers={"Content-Type": "application/json"},
            header_auth=True,
        )

    elif args.request_method == "DELETE":
        r = session.delete(url)
    if r.status_code == 200 or r.status_code == 201:
        if re.search("json", r.headers["content-type"], flags=0):
            print(r.json())
        else:
            print(r.text)
    elif r.status_code == 400:
        print("Bad request")
        print(r.json()["message"])
    elif r.status_code == 401:
        if re.search("unauthorized", r.json()["message"]):
            print("Access denied - client is unauthorized")
            return
        else:
            if re.search("verification", r.json()["message"]):
                print(r.json())
                return
            print("Invalid session token, requesting new one...\n")
            (token, secret) = get_session_token(None, None)
            get_route(route, token, secret, args)
    else:
        print("Error:")
        print(r.text)
    return


def get_request_token(service):
    print("Getting request token...")

    r = service.get_raw_request_token(params={"oauth_callback": "oob"})
    if r.status_code == 200:
        print("Success:")
        token = r.json()["oauth_token"]
        secret = r.json()["oauth_token_secret"]
        print("Request Token:        " + token)
        print("Request Token Secret: " + secret + "\n")
        write_token("request_token", token, secret)
        return (token, secret)
    else:
        print("Failed:")
        print(r.json()["message"])
    return


def get_access_token(request_token, request_secret, service):
    if os.path.exists(ACCESS_TOKEN_PATH):
        os.remove(ACCESS_TOKEN_PATH)
    if not request_token or not request_secret:
        (request_token, request_secret) = get_request_token(service)
    print("Now log in at " + AUTHORIZE_URL + "&oauth_token=" + request_token + "\n")
    verifier = input("Please enter verification code: ")
    r = service.get_raw_access_token(
        request_token, request_secret, params={"oauth_verifier": verifier}
    )
    if r.status_code == 200:
        print("Success:")
        token = r.json()["oauth_token"]
        secret = r.json()["oauth_token_secret"]
        print("Access Token:        " + token)
        print("Access Token Secret: " + secret + "\n")
        print("This access token will not expire but may be revoked")
        print("by the user or the service provider. It may be used to")
        print("obtain temporary session tokens.")
        write_token(ACCESS_TOKEN_PATH, token, secret)
        return (token, secret)
    else:
        print("Failed:")
        print(r.json()["message"])
        os._exit(1)
    return


def get_session_token(access_token, access_secret):
    if os.path.exists(SESSION_TOKEN):
        os.remove(SESSION_TOKEN)
    if not access_token or not access_secret:
        (access_token, access_secret) = retrieve_token(ACCESS_TOKEN_PATH)
        if not access_token or not access_secret:
            raise Exception(f'Access token {ACCESS_TOKEN_PATH} is missing - please check.')
            # (access_token, access_secret) = get_access_token(None, None)

    print("Now requesting session token using access token...\n")
    session_request = OAuth1Session(
        CONSUMER_KEY,
        CONSUMER_SECRET,
        access_token=access_token,
        access_token_secret=access_secret,
    )
    url = SESSION_REST_URL + "/oauth/get_session_token"

    r = session_request.get(url)
    if r.status_code == 200:
        session_token = r.json()["oauth_token"]
        session_secret = r.json()["oauth_token_secret"]
        print("Session Token:        " + session_token)
        print("Session Token Secret: " + session_secret + "\n")
        print("This session token will expire in 12 hours (default).")
        print("It should be used with the secret to sign any requests")
        print("to the API.")
        write_token(SESSION_TOKEN, session_token, session_secret)
        return (session_token, session_secret)
    else:
        print("Failed:")
        print(r.json()["message"])
        os._exit(1)
    return


# def prompt():
#     if not args.prompt:
#         return
#     input("Press any key to continue...")
#     return


def write_token(token_type, token, secret):
    token_file = open(token_type, "w")
    token_file.write("token=" + token + "\n")
    token_file.write("secret=" + secret + "\n")
    token_file.close()


if __name__ == "__main__":
    main()
