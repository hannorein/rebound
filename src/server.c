/**
 * @file    server.c
 * @brief   Opens a webserver to allow for platform independent visualization.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details These functions provide real time visualizations
 * using OpenGL. Part of the code is by Dave O'Hallaron, Carnegie Mellon (tiny.c).
 * 
 * @section LICENSE
 * Copyright (c) 2023 Hanno Rein, Dave O'Hallaron, Carnegie Mellon
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef SERVER

#include <stdio.h>
#ifdef _MSC_VER 
//not #if defined(_WIN32) || defined(_WIN64) because we have strncasecmp in mingw
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif
#ifdef _WIN32
#include <WS2tcpip.h>
#include <tchar.h>
#include <io.h>
#define F_OK 0
#define access _access
#pragma comment(lib, "ws2_32.lib")
#else // _WIN32
#include <unistd.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/mman.h>
#include <netinet/in.h>
#endif // _WIN32
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include "rebound.h"

#define BUFSIZE 1024

const char* reb_server_header =
        "HTTP/1.1 200 OK\n"
        "Server: REBOUND Webserver\n"
      //"Access-Control-Allow-Origin: *\n"
      //"Cross-Origin-Opener-Policy: cross-origin\n"
        "Content-type: text/html\n"
        "\r\n";
const char* reb_server_header_png =
        "HTTP/1.1 200 OK\n"
        "Server: REBOUND Webserver\n"
        "Content-type: image/png\n"
        "\r\n";

#ifdef _WIN32
int sendBytes(SOCKET s, const void * buffer, int buflen){
    int total = 0;
    char *pbuf = (char*) buffer;
    while (buflen > 0) {
        int iResult = send(s, pbuf, buflen, 0);
        if (iResult < 0) {
            if (WSAGetLastError() == WSAEWOULDBLOCK) {
                // optionally use select() to wait for the
                // socket to have more space to write before
                // calling send() again...
                continue;
            }

            printf("send error: %d\n", WSAGetLastError());
            return SOCKET_ERROR;
        } else if (iResult == 0) {
            printf("disconnected\n");
            return 0;
        } else {
            pbuf += iResult;
            buflen -= iResult;
            total += iResult;
        }
    }

    return total;
}
#endif // _WIN32


#ifdef _WIN32
static void reb_server_cerror(SOCKET clientS, char cause[]){
#else //_WIN32
static void reb_server_cerror(FILE *stream, char *cause){
#endif //_WIN32
    char* buf = NULL;
    asprintf(&buf,  "HTTP/1.1 501 Not Implemented\n"
                    "Content-type: text/html\n"
                    "\n"
                    "<html><title>REBOUND Webserver Error</title>"
                    "<body>\n"
                    "<h1>Error</h1>\n"
                    "<p>%s</p>\n"
                    "<hr><em>REBOUND Webserver</em>\n"
                    "</body></html>\n"
                    , cause);
    printf("\nREBOUND Webserver error: %s\n", cause);
#ifdef _WIN32
    sendBytes(clientS, buf, strlen(buf));
    closesocket(clientS); // close the Client Socket now that our Work is Complete.
#else //_WIN32
    fwrite(buf, 1, strlen(buf), stream);
#endif //_WIN32
    free(buf);
}


void* reb_server_start(void* args){
    struct reb_server_data* data = (struct reb_server_data*)args;
    struct reb_simulation* r = data->r;

    if (access("rebound.html", F_OK)) {
        reb_simulation_warning(r, "File rebound.html not found in current directory. Attempting to download it from github.");
        char curl_cmd[] = "curl -s --output rebound.html https://gist.githubusercontent.com/hannorein/8ad8abc9bb5d2a76cf8da2b104ef8c79/raw/ce41328f2121b9a4d4578cb592684f4c276ba9ed/rebound.html";
        system(curl_cmd);
        if (access("rebound.html", F_OK)) {
            reb_simulation_warning(r, "Automatic download failed. Manually download the file from github and place it in the current directory to enable browser based visualization.");
        }else{
            printf("Success: rebound.html downloaded.\n");
        }
    }

#ifndef _WIN32
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);


    /* variables for connection management */
    int childfd;           /* child socket */
    unsigned int clientlen;         /* byte size of client's address */
    struct hostent *hostp; /* client host info */
    int optval;            /* flag value for setsockopt */
    struct sockaddr_in serveraddr; /* server's addr */
    struct sockaddr_in clientaddr; /* client addr */

    /* variables for connection I/O */
    FILE *stream;          /* stream version of childfd */
    char buf[BUFSIZE];     /* message buffer */
    char method[BUFSIZE];  /* request method */
    char uri[BUFSIZE];     /* request uri */
    char version[BUFSIZE]; /* request method */

    /* open socket descriptor */
    data->socket = socket(AF_INET, SOCK_STREAM, 0);
    if (data->socket < 0)
        reb_exit("ERROR opening socket");

    /* allows us to restart server immediately */
    optval = 1;
    setsockopt(data->socket, SOL_SOCKET, SO_REUSEADDR,
            (const void *)&optval , sizeof(int));

    /* bind port to socket */
    memset((char *) &serveraddr, 0, sizeof(serveraddr));
    serveraddr.sin_family = AF_INET;
    serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
    serveraddr.sin_port = htons((unsigned short)data->port);
    if (bind(data->socket, (struct sockaddr *) &serveraddr, sizeof(serveraddr)) < 0){
        char error_msg[BUFSIZE];
        snprintf(error_msg, BUFSIZE, "Error binding to port %d. Port might be in use.\n", data->port);
        reb_simulation_error(r, error_msg);
        data->ready = -1;
        return PTHREAD_CANCELED;
    }

    /* get us ready to accept connection requests */
    if (listen(data->socket, 5) < 0) /* allow 5 requests to queue up */
        reb_exit("ERROR on listen");

    printf("REBOUND Webserver listening on http://localhost:%d ...\n",data->port);
    /*
     * main loop: wait for a connection request, parse HTTP,
     * serve requested content, close connection.
     */
    clientlen = sizeof(clientaddr);
    while (1) {
        /* wait for a connection request */
        data->ready = 1;
        childfd = accept(data->socket, (struct sockaddr *) &clientaddr, &clientlen);
        if (childfd < 0) { // Accept will fail if main thread is closing socket.
            return 1;
        }

        /* determine who sent the message */
        hostp = gethostbyaddr((const char *)&clientaddr.sin_addr.s_addr,
                sizeof(clientaddr.sin_addr.s_addr), AF_INET);
        if (hostp == NULL)
            reb_exit("ERROR on gethostbyaddr");

        /* open the child socket descriptor as a stream */
        if ((stream = fdopen(childfd, "r+")) == NULL)
            reb_exit("ERROR on fdopen");

        /* get the HTTP request line */
        char* request = fgets(buf, BUFSIZE, stream);
        if (!request){
            reb_server_cerror(stream, "Did not get request.");
            fclose(stream);
            close(childfd);
            continue;
        }
        sscanf(buf, "%s %s %s\n", method, uri, version);

        /* only support the GET method */
        if (strcasecmp(method, "GET")) {
            reb_server_cerror(stream, "Only GET is implemented.");
            fclose(stream);
            close(childfd);
            continue;
        }
           
        /* read (and ignore) the HTTP headers */
        fgets(buf, BUFSIZE, stream);
        while(strcmp(buf, "\r\n")) {
            fgets(buf, BUFSIZE, stream);
        }

        if (!strcasecmp(uri, "/simulation")) {
            char* bufp = NULL;
            size_t sizep;
            data->need_copy = 1;
            pthread_mutex_lock(&(data->mutex));
            reb_simulation_save_to_stream(r, &bufp,&sizep);
            data->need_copy = 0;
            pthread_mutex_unlock(&(data->mutex));
            fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
            fwrite(bufp, 1, sizep, stream);
            free(bufp);
        }else if (!strncasecmp(uri, "/keyboard/",10)) {
            int key = 0;
            sscanf(uri, "/keyboard/%d", &key);
            switch (key){
                case 'Q':
                    data->r->status = REB_STATUS_USER;
                    fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
                    fprintf(stream, "ok.\n");
                    break;
                case ' ':
                    if (data->r->status == REB_STATUS_PAUSED){
                        printf("Resume.\n");
                        data->r->status = REB_STATUS_RUNNING;
                    }else{
                        printf("Pause.\n");
                        data->r->status = REB_STATUS_PAUSED;
                    }
                    fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
                    fprintf(stream, "ok.\n");
                    break;
                default:
                    reb_server_cerror(stream, "Unsupported key received.");
                    break;
            }
        }else if (!strcasecmp(uri, "/") || !strcasecmp(uri, "/index.html") || !strcasecmp(uri, "/rebound.html")) {
            struct stat sbuf;
            if (stat("rebound.html", &sbuf) < 0) {
                reb_server_cerror(stream, "rebound.html not found in current directory. Try `make rebuund.html`.");
            }else{
                fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
                int fd = open("rebound.html", O_RDONLY);
                void* p = mmap(0, sbuf.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
                fwrite(p, 1, sbuf.st_size, stream);
                munmap(p, sbuf.st_size);
            }
        }else if (!strcasecmp(uri, "/favicon.ico")) {
                fwrite(reb_server_header_png, 1, strlen(reb_server_header_png), stream);
                fwrite(reb_favicon_png,1, reb_favicon_len, stream);
        }else{
            reb_server_cerror(stream, "Unsupported URI.");
        }

        /* clean up */
        fflush(stream);
        fclose(stream);
        close(childfd);

    }
    printf("Server shutting down...\n");
    return PTHREAD_CANCELED;

#else // _WIN32


    WSADATA wsa;
    struct sockaddr_in server;
    SOCKET clientS;
    char request[BUFSIZE];
    char method[BUFSIZE];
    char uri[BUFSIZE];
    char version[BUFSIZE];

    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) {
        printf("Winsock startup failed");
        exit(1);
    }

    data->socket = socket(AF_INET, SOCK_STREAM, 0);
    if (data->socket == INVALID_SOCKET) {
        printf("Socket error\n");
        exit(1);
    }

    server.sin_family = AF_INET;
    server.sin_port = htons(data->port);
    InetPton(AF_INET, _T("0.0.0.0"), &server.sin_addr);

    int ret_bind = bind(data->socket, (struct sockaddr*)&server, sizeof(server)); // binding the Host Address and Port Number
    if (ret_bind) {
        char error_msg[BUFSIZE];
        snprintf(error_msg, BUFSIZE, "Error binding to port %d. Port might be in use.\n", data->port);
        reb_simulation_error(r, error_msg);
        data->ready = -1;
        return 1;
    }

    int ret_listen = listen(data->socket, AF_INET);
    if (ret_listen){
        printf("Listen error\n");
        exit(1);
    }
    
    printf("REBOUND Webserver listening on http://localhost:%d ...\n",data->port);

    while(1){
        data->ready = 1;
        clientS = accept(data->socket, NULL, NULL);
        if (clientS == INVALID_SOCKET) { // Accept will fail if main thread is closing socket.
            return 1;
        } 
        // Receive request. Ideally we should check for new line and read more bytes if needed.
        recv(clientS, request, BUFSIZE, 0);

        sscanf(request, "%s %s %s\n", method, uri, version);

        /* only support the GET method */
        if (strcasecmp(method, "GET")) {
            reb_server_cerror(clientS, "Method not Implemented");
            continue;
        }
        
        if (!strcasecmp(uri, "/simulation")) {
            char* bufp = NULL;
            size_t sizep;
            data->need_copy = 1;
            //pthread_mutex_lock(&(data->mutex)); TODO!!
            reb_simulation_save_to_stream(r, &bufp,&sizep);
            data->need_copy = 0;
            //pthread_mutex_unlock(&(data->mutex));
            sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
            sendBytes(clientS, bufp, sizep);
            free(bufp);
        }else if (!strncasecmp(uri, "/keyboard/",10)) {
            int key = 0;
            sscanf(uri, "/keyboard/%d", &key);
            switch (key){
                case 'Q':
                    data->r->status = REB_STATUS_USER;
                    break;
                case ' ':
                    if (data->r->status == REB_STATUS_PAUSED){
                        printf("Resume.\n");
                        data->r->status = REB_STATUS_RUNNING;
                    }else{
                        printf("Pause.\n");
                        data->r->status = REB_STATUS_PAUSED;
                    }
                    break;
                default:
                    reb_server_cerror(clientS, "Unknown key received.");
                    continue;
                    break;
            }
        }else if (!strcasecmp(uri, "/") || !strcasecmp(uri, "/index.html") || !strcasecmp(uri, "/rebound.html")) {
            FILE *f = fopen("rebound.html", "rb");
            if (f){
                fseek(f, 0, SEEK_END);
                long fsize = ftell(f);
                fseek(f, 0, SEEK_SET);
                char *buf = malloc(fsize);
                fread(buf, fsize, 1, f);
                fclose(f);
                sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
                sendBytes(clientS, buf, fsize); 
                free(buf);
            }else{
                reb_server_cerror(clientS, "rebound.html not found in current directory. Try `make rebuund.html`.");
                continue;
            }
        }else if (!strcasecmp(uri, "/favicon.ico")) {
                sendBytes(clientS, reb_server_header_png, strlen(reb_server_header_png)); 
                sendBytes(clientS, reb_favicon_png, reb_favicon_len); 
        }else{
            reb_server_cerror(clientS, "Unsupported request.");
            printf("URI: %s\n",uri);
            continue;
        }

        closesocket(clientS);
    }
    WSACleanup();
    return NULL;
#endif // _WIN32
}


#endif // SERVER
