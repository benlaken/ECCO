#!/usr/bin/env bash
case "$1" in
    dev)
        docker-compose -f docker-compose-process.yml build && docker-compose -f docker-compose-process.yml up
        ;;
    test)
        docker run --name postgis -e POSTGRES_PASSWORD=password123 -d mdillon/postgis
        docker run -it --link postgis:postgres --rm postgres \
          sh -c 'exec psql -h "$POSTGRES_PORT_5432_TCP_ADDR" -p "$POSTGRES_PORT_5432_TCP_PORT" -U postgres'
        ;;
    production)
        echo "This does nothing yet"
        ;;
    *)
        exec "$@"
esac
