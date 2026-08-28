function connection = timeout_connection(connection)
%TIMEOUT_CONNECTION Mark a connection as stopped by a declared runtime cap.
connection.success = false;
connection.totalIterations = Inf;
connection.tauResolved = Inf;
connection.diagnostics.runtimeLimitExceeded = true;
connection.diagnostics.failureCategory = "runtime_limit_exceeded";
end
