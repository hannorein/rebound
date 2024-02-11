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

#include "rebound.h"

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
#include <errno.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/mman.h>
#include <netinet/in.h>
#endif // _WIN32
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>



#define BUFSIZE 1024

const char* reb_server_header =
        "HTTP/1.1 200 OK\n"
        "Server: REBOUND Webserver\n"
        "Cache-Control: no-cache, no-store, must-revalidate\n"
        "Pragma: no-cache\n"
        "Expires: 0\n"
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

static const unsigned char base64_table[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/**
 * base64_decode - Base64 decode
 * @src: Data to be decoded
 * @len: Length of the data to be decoded
 * @out_len: Pointer to output length variable
 * Returns: Allocated buffer of out_len bytes of decoded data,
 * or %NULL on failure
 *
 * Caller is responsible for freeing the returned buffer.
 *
 * Source: https://web.mit.edu/freebsd/head/contrib/wpa/src/utils/base64.c
 */
static unsigned char * base64_decode(const unsigned char *src, size_t len, size_t *out_len) {
	unsigned char dtable[256], *out, *pos, block[4], tmp;
	size_t i, count, olen;
	int pad = 0;

	memset(dtable, 0x80, 256);
	for (i = 0; i < sizeof(base64_table) - 1; i++)
		dtable[base64_table[i]] = (unsigned char) i;
	dtable['='] = 0;

	count = 0;
	for (i = 0; i < len; i++) {
		if (dtable[src[i]] != 0x80)
			count++;
	}

	if (count == 0 || count % 4)
		return NULL;

	olen = count / 4 * 3;
	pos = out = malloc(olen);
	if (out == NULL)
		return NULL;

	count = 0;
	for (i = 0; i < len; i++) {
		tmp = dtable[src[i]];
		if (tmp == 0x80)
			continue;

		if (src[i] == '=')
			pad++;
		block[count] = tmp;
		count++;
		if (count == 4) {
			*pos++ = (block[0] << 2) | (block[1] >> 4);
			*pos++ = (block[1] << 4) | (block[2] >> 2);
			*pos++ = (block[2] << 6) | block[3];
			count = 0;
			if (pad) {
				if (pad == 1)
					pos--;
				else if (pad == 2)
					pos -= 2;
				else {
					/* Invalid padding */
					free(out);
					return NULL;
				}
				break;
			}
		}
	}

	*out_len = pos - out;
	return out;
}


#ifndef _WIN32
void* reb_server_start(void* args){
#else //_WIN32
DWORD WINAPI reb_server_start(void* args){
#endif // _WIN32
    struct reb_server_data* data = (struct reb_server_data*)args;
    struct reb_simulation* r = data->r;

    if (access("rebound.html", F_OK)) {
        reb_simulation_warning(r, "File rebound.html not found in current directory. Attempting to download it from github.");
        char curl_cmd[] = "curl -L -s --output rebound.html https://github.com/hannorein/rebound/releases/latest/download/rebound.html";
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
    int childfd;                   /* child socket */
    unsigned int clientlen;        /* byte size of client's address */
    int optval;                    /* flag value for setsockopt */
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
    serveraddr.sin_port = htons(data->port);
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
            return PTHREAD_CANCELED;
        }

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
        if (strcasecmp(method, "GET") && strcasecmp(method, "POST")) {
            reb_server_cerror(stream, "Only GET+POST are implemented.");
            fclose(stream);
            close(childfd);
            continue;
        }
           
        /* read (and ignore) the HTTP headers */
        fgets(buf, BUFSIZE, stream);
        unsigned long content_length = 0;
        while(strcmp(buf, "\r\n")) {
            char cl[BUFSIZE];
            int ni = sscanf(buf, "Content-Length: %s\n", cl);
            if (ni){
                content_length = strtol(cl,NULL,10);
            }
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
            data->need_copy = 1;
            pthread_mutex_lock(&(data->mutex));
            int skip_default_keys = 0;
            if (r->key_callback){
                skip_default_keys = r->key_callback(r, key);
            } 
            data->need_copy = 0;
            pthread_mutex_unlock(&(data->mutex));
            if (!skip_default_keys){
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
                        }else if (data->r->status == REB_STATUS_RUNNING){
                            printf("Pause.\n");
                            data->r->status = REB_STATUS_PAUSED;
                        }
                        fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
                        fprintf(stream, "ok.\n");
                        break;
                    case 264: // arrow down
                        if (data->r->status == REB_STATUS_PAUSED){
                            data->r->status = REB_STATUS_SINGLE_STEP;
                            printf("Step.\n");
                        }
                        fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
                        fprintf(stream, "ok.\n");
                        break;
                    case 267: // page down
                        if (data->r->status == REB_STATUS_PAUSED){
                            data->r->status = REB_STATUS_SINGLE_STEP - 50;
                            printf("50 steps.\n");
                        }
                        fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
                        fprintf(stream, "ok.\n");
                        break;
                    default:
                        // reb_server_cerror(stream, "Unsupported key received.");
                        break;
                }
            }else{
                fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
                fprintf(stream, "ok.\n");
            }
        }else if (!strcasecmp(uri, "/") || !strcasecmp(uri, "/index.html") || !strcasecmp(uri, "/rebound.html")) {
            struct stat sbuf;
            if (stat("rebound.html", &sbuf) < 0) {
                reb_server_cerror(stream, "rebound.html not found in current directory. Try `make rebound.html`.");
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
        }else if (!strcasecmp(uri, "/screenshot")) {
            data->need_copy = 1;
            pthread_mutex_lock(&(data->mutex));
            
            if (content_length==0){
                printf("Received screenshot with size zero.");
                goto screenshot_finish;
            }
            if (r->status != REB_STATUS_SCREENSHOT){
                printf("Received screenshot but did not expect one.\n");
                goto screenshot_finish;
            }
            if (data->screenshot) {
                printf("Unable to receive screenshot as previous screenshot not freed.\n");
                goto screenshot_finish;
            }

            char* dataURL = malloc(content_length);
            int rc = fread(dataURL, content_length, 1,  stream);

            if (rc!=1){
                printf("Error while reading screenshot data.\n");
                free(dataURL);
                goto screenshot_finish;
            }

            int rc_len = strlen(dataURL)+1;
            char* base64 = strchr(dataURL, ',');
            if (content_length != rc_len){
                printf("Received screenshot with incorrect size.\n");
                free(dataURL);
                goto screenshot_finish;
            }
            if (!base64){
                printf("Unable to decode received screenshot. Data not in dataURL format.\n");
                free(dataURL);
                goto screenshot_finish;
            }
            data->screenshot = base64_decode((unsigned char*)base64+1, strlen(base64+1), &data->N_screenshot);
            if (!data->screenshot){
                printf("An error occured while decoding the screenshot.\n");
            }
            data->r->status = REB_STATUS_PAUSED;
            free(dataURL);
screenshot_finish:
            data->need_copy = 0;
            pthread_mutex_unlock(&(data->mutex));
            fwrite(reb_server_header, 1, strlen(reb_server_header), stream);
            fprintf(stream, "ok.\n");
        }else{
            reb_server_cerror(stream, "Unsupported URI.");
            printf("URI: %s\n", uri);
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
        // Receive entire request.
        int recN = 0;
        char* recbuf = malloc(BUFSIZE);
        int recbufN = 0;
        while((recN = recv(clientS, recbuf+recbufN, BUFSIZE, 0))>0){
            recbufN += recN;
            if (recN<BUFSIZE) break;
            recbuf = realloc(recbuf, recbufN+BUFSIZE);
        }

        // Get method and uri
        sscanf(recbuf, "%s %s %s\n", method, uri, version);
        if (strcasecmp(method, "GET") && strcasecmp(method, "POST")) {
            reb_server_cerror(clientS, "Method not Implemented");
            continue;
        }
        
        /* read (and ignore) the HTTP headers */
        char* curLine = recbuf;
        unsigned long content_length = 0;
        while(curLine){
            char* nextLine = strchr(curLine, '\n');
            if (nextLine) *nextLine = '\0';
            char cl[BUFSIZE];
            int ni = sscanf(curLine, "Content-Length: %s\n", cl);
            if (ni){
                content_length = strtol(cl,NULL,10);
            }
            if (nextLine){
                *nextLine = '\n';
                curLine = nextLine+1;
            }else{
                break;
            }
        }
        
        // Only post data is needed, otherwise free here
        if (strcasecmp(method, "POST")) {
            free(recbuf);
        }
        
        if (!strcasecmp(uri, "/simulation")) {
            char* bufp = NULL;
            size_t sizep;
            data->need_copy = 1;
            WaitForSingleObject(data->mutex, INFINITE);
            reb_simulation_save_to_stream(r, &bufp,&sizep);
            data->need_copy = 0;
            ReleaseMutex(data->mutex);
            sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
            sendBytes(clientS, bufp, sizep);
            free(bufp);
        }else if (!strncasecmp(uri, "/keyboard/",10)) {
            int key = 0;
            const char* ok = "ok.";
            sscanf(uri, "/keyboard/%d", &key);
            int skip_default_keys = 0;
            data->need_copy = 1;
            WaitForSingleObject(data->mutex, INFINITE);
            if (r->key_callback){
                skip_default_keys = r->key_callback(r, key);
            } 
            data->need_copy = 0;
            ReleaseMutex(data->mutex);
            if (!skip_default_keys){
                switch (key){
                    case 'Q':
                        data->r->status = REB_STATUS_USER;
                        sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
                        sendBytes(clientS, ok, strlen(ok));
                        break;
                    case ' ':
                        if (data->r->status == REB_STATUS_PAUSED){
                            printf("Resume.\n");
                            data->r->status = REB_STATUS_RUNNING;
                        }else if (data->r->status == REB_STATUS_RUNNING){
                            printf("Pause.\n");
                            data->r->status = REB_STATUS_PAUSED;
                        }
                        sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
                        sendBytes(clientS, ok, strlen(ok));
                        break;
                    case 264: // down arrow
                        if (data->r->status == REB_STATUS_PAUSED){
                            data->r->status = REB_STATUS_SINGLE_STEP;
                            printf("Step.\n");
                        }
                        sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
                        sendBytes(clientS, ok, strlen(ok));
                        break;
                    case 267: // page down
                        if (data->r->status == REB_STATUS_PAUSED){
                            data->r->status = REB_STATUS_SINGLE_STEP - 50;
                            printf("50 step.\n");
                        }
                        sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
                        sendBytes(clientS, ok, strlen(ok));
                        break;
                    default:
                        // reb_server_cerror(clientS, "Unknown key received.");
                        continue;
                        break;
                }
            }else{
                sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
                sendBytes(clientS, ok, strlen(ok));
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
                reb_server_cerror(clientS, "rebound.html not found in current directory. Try `make rebound.html`.");
                continue;
            }
        }else if (!strcasecmp(uri, "/favicon.ico")) {
                sendBytes(clientS, reb_server_header_png, strlen(reb_server_header_png)); 
                sendBytes(clientS, reb_favicon_png, reb_favicon_len); 
        }else if (!strcasecmp(uri, "/screenshot")) {
            data->need_copy = 1;
            WaitForSingleObject(data->mutex, INFINITE);
            if (content_length==0){
                printf("Received screenshot with size zero.");
                goto screenshot_finish;
            }
            if (r->status != REB_STATUS_SCREENSHOT){
                printf("Received screenshot but did not expect one.\n");
                goto screenshot_finish;
            }
            if (data->screenshot) {
                printf("Unable to receive screenshot as previous screenshot not freed.\n");
                goto screenshot_finish;
            }

            char* dataURL = curLine; // Memory!

            int rc_len = strlen(dataURL)+1;
            char* base64 = strchr(dataURL, ',');
            if (content_length != rc_len){
                printf("Received screenshot with incorrect size.\n");
                goto screenshot_finish;
            }
            if (!base64){
                printf("Unable to decode received screenshot. Data not in dataURL format.\n");
                goto screenshot_finish;
            }
            data->screenshot = base64_decode((unsigned char*)base64+1, strlen(base64+1), &data->N_screenshot);
            if (!data->screenshot){
                printf("An error occured while decoding the screenshot.\n");
            }
            data->r->status = REB_STATUS_PAUSED;
screenshot_finish:
            free(recbuf);
            data->need_copy = 0;
            ReleaseMutex(data->mutex);
            const char* ok = "ok.";
            sendBytes(clientS, reb_server_header, strlen(reb_server_header)); 
            sendBytes(clientS, ok, strlen(ok));
        }else{
            reb_server_cerror(clientS, "Unsupported request.");
            printf("URI: %s\n",uri);
            continue;
        }

        closesocket(clientS);
    }
    WSACleanup();
    return 0;
#endif // _WIN32
}

#endif // SERVER


int reb_simulation_start_server(struct reb_simulation* r, int port){    
#ifdef SERVER
    if (port){
        if (r->server_data){
            reb_simulation_error(r,"Server already started.");
            return -1;
        }
        r->server_data = calloc(sizeof(struct reb_server_data),1);
        r->server_data->r = r;
        r->server_data->port = port;
#ifdef _WIN32
        r->server_data->mutex = CreateMutex(NULL, FALSE, NULL);
        HANDLE thread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)reb_server_start, r->server_data, 0, NULL);
#else // _WIN32
        if (pthread_mutex_init(&(r->server_data->mutex), NULL)){
            reb_simulation_error(r,"Mutex creation failed.");
            return -1;
        }
        int ret_create = pthread_create(&(r->server_data->server_thread),NULL,reb_server_start,r->server_data);
        if (ret_create){
            reb_simulation_error(r, "Error creating server thread.");
            return -1;
        }
#endif // _WIN32
        int maxwait = 100;
        while (r->server_data->ready==0 && maxwait){
            usleep(10000);
            maxwait--;
        }
        if (r->server_data->ready==0){
            reb_simulation_warning(r, "Server did not start immediately. This might just take a little bit longer.");
        }
        return 0;
    }else{
        reb_simulation_error(r, "Cannot start server. Invalid port.");
        return -1;
    }
#else // SERVER
#ifndef SERVERHIDEWARNING
    reb_simulation_error(r, "REBOUND has been compiled without SERVER support.");
#endif // SERVERHIDEWARNING
    return -1;
#endif // SERVER
}

void reb_simulation_stop_server(struct reb_simulation* r){    
#ifdef SERVER
    if (r==NULL) return;
    if (r->server_data){
#ifdef _WIN32
        closesocket(r->server_data->socket); // Will cause thread to exit.
#else // _WIN32
        close(r->server_data->socket); // Will cause thread to exit.
        int ret_cancel = pthread_cancel(r->server_data->server_thread);
        if (ret_cancel==ESRCH){
            printf("Did not find server thread while trying to cancel it.\n");
        }
        void* retval = 0;
        pthread_join(r->server_data->server_thread, &retval);
        if (retval!=PTHREAD_CANCELED){
            printf("An error occured while cancelling server thread.\n");
        }
#endif // _WIN32
        free(r->server_data);
        r->server_data = NULL;
    }
#endif //SERVER
}



